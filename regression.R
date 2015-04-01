
library(optparse)
source('~/softwares/my/ml_util.R')

opt <- interface()

if (opt$debug) save.image('debug.Rdata')

if (opt$verbose) {
    cat("Running command with args:\n",
        paste(commandArgs(), collapse = " "),
        '\n')
}


if (opt$split <= 0 & opt$split > 1) {
    stop("The split arg should be greater than 0 and not greater than 1.")
}

if(is.null(opt$models)) {
    models <- regression
} else {
    models <- strsplit(opt$models, ',')[[1]]
}

library(caret)
if (opt$cores > 1) {
    library(doMC)
    registerDoMC(opt$cores)
}

meta <- read.table.x(opt$metadata, quote='')

meta.col <- colnames(meta)


if (is.null(opt$fields)) {
    ## stop("No field is provided to do regression on.")
    outcome.col <- names(meta)
    boring <- c("#SampleID",
                "BarcodeSequence",
                "LinkerPrimerSequence",
                "TARGET_SUBFRAGMENT",
                "ASSIGNED_FROM_GEO",
                "EXPERIMENT_CENTER",
                "RUN_PREFIX",
                "TAXON_ID",
                "ILLUMINA_TECHNOLOGY",
                "COMMON_NAME",
                "EXTRACTED_DNA_AVAIL_NOW",
                "SAMPLE_CENTER",
                "STUDY_CENTER",
                "Description")
    outcome.col <- outcome.col[! outcome.col %in% boring]
} else {
    outcome.col <- strsplit(opt$fields, ',')[[1]]
    x <- which(! outcome.col %in% meta.col)
    if (length(x) > 0) {
        stop("Field(s) ", paste(outcome.col[x], collapse=','), " do not exist in meta data")
    }
}

## extract part of the samples by their meta data
if (! is.null(opt$category)) {
    ## e.g. "SITE::nostril,skin:_:SEX::male"
    extract <- strsplit(opt$category, ':_:')[[1]]
    extract <- strsplit(extract, '::')
    for (x in extract) {
        if (! x[1] %in% meta.col)
            stop("The field ", x[1], " does not exist in meta data")
        i <- meta[[ x[1] ]]
        j <- strsplit(x[2], ',')[[1]]
        if (! all(j %in% i)) {
            ## insanity check to avoid typos
            stop("You specified non-existing values for field ", x[1], " in meta data")
        }
        meta <- meta[ i %in% j, ]
    }
}


if (! is.null(opt$numeric)) {
    ## e.g. "PH::6,12:_:TEMP::,32"
    extract <- strsplit(opt$numeric, ':_:')[[1]]
    extract <- strsplit(extract, '::')
    for (x in extract) {
        if (! x[1] %in% meta.col)
            stop("The field ", x[1], " does not exist in meta data")
        ## in case there is None, NA, etc in the column (R will read it into character
        ## instead of numerical)
        i <- as.numeric(as.character(meta[[ x[1] ]]))
        j <- as.numeric(strsplit(x[2], ',')[[1]])
        ## NA & TRUE -> NA
        ## NA & FALSE -> FALSE
        n <- rep(T, nrow(meta))
        if (! is.na(j[1]))
            n <- n & i > j[1]
        if (! is.na(j[2]))
            n <- n & i < j[2]
        meta <- meta[n, ]
    }
}


otus <- read.table.x(opt$input_otu_table)


tax.16s <- otus[, length(otus)]
tax.16s <- gsub("^Root; ", "", tax.16s)
## insert a newline for every three levels of taxonomy
tax.16s <- gsub("([^;]*); ([^;]*); ([^;]*); ", '\\1; \\2; \\3\n', tax.16s)
names(tax.16s) <- otus[[1]]


## remove the 6-digit suffix of the sample IDs in the mapping file.
## meta.sid <- gsub(".[0-9]{6}$", "", as.character(meta[[1]]))
meta.sid <- as.character(meta[[1]])
rownames(meta) <- meta.sid
sample.ids <- intersect(meta.sid, colnames(otus))
meta <- meta[sample.ids, ]


rownames(otus) <- otus[[1]]
otus <- data.frame(t(otus[, sample.ids]), check.names=FALSE)


## add the categorical fields in the meta data as predictors
if (! is.null(opt$add_numeric)) {
    add.pred <- strsplit(opt$add_numeric, ',', fixed=TRUE)[[1]]
    not.in <- which(! add.pred %in% colnames(meta))
    if (length(not.in) > 0) {
        stop(paste(c(add.pred[not.in], "not in the meta data!!!"), collapse=' '))
    }
    otus <- cbind(as.numeric(meta[, add.pred]), otus)
    colnames(otus)[1:length(add.pred)] <- add.pred
}
## add numeric fields as predictors
if (! is.null(opt$add_category)) {
    add.pred <- strsplit(opt$add_category, ',', fixed=TRUE)[[1]]
    not.in <- which(! add.pred %in% colnames(meta))
    if (length(not.in) > 0) {
        stop(paste(c(add.pred[not.in], "not in the meta data!!!"), collapse=' '))
    }
    if (opt$debug) save.image('debug.Rdata')
    x <- meta[, add.pred, drop=FALSE]
    dummy <- dummyVars(~., data=x)
    otus <- cbind(predict(dummy, x), otus)
}

pdf(sprintf("%s.pdf", opt$output))
save.image(sprintf("%s.Rdata", opt$output))

min.sample.size <- 12
accuracies <- data.frame()
top.features <- data.frame()
big.tuned.list <- list()
for(label in outcome.col) {
    if(opt$verbose)
        cat("=========", label, ":\n")

    outcome <- as.numeric(as.character(meta[[label]]))

    if(opt$verbose) {
        cat("---- a glimpse of outcome:\n")
        print(head(outcome, n=30))
        print(summary(outcome))
    }

    ## remove NA values
    outcome.na <- is.na(outcome)
    ## if more than half of the samples are not numeric
    if(sum(outcome.na) > 0.5 * length(outcome)) {
        if(opt$verbose)
            cat("outcome has less than half of numeric values. skip it.\n")
        next
    }

    outcome <- outcome[! outcome.na]
    if (length(outcome) < min.sample.size)
        stop("There should be more than ", min.sample.size, " samples.")

    ## if there less than 3 uniq values in this category
    if(length(unique(outcome)) < 2) {
        if(opt$verbose)
            cat("outcome has less than 2 distinctive values. skip it.\n")
        next
    }

    train.set <- otus[! outcome.na, ]


    if (opt$split < 1) {
        set.seed(1)
        training.rows <- createDataPartition(outcome, p=opt$split, list=F)
    } else {
        training.rows <- 1:length(outcome)
    }

    train.full <- train.set[training.rows, ]
    test.full <- train.set[-training.rows, ]
    train.outcome <- outcome[training.rows]
    test.outcome <- outcome[-training.rows]

    nzv <- nearZeroVar(train.full)
    if (length(nzv) > 0) {
        train.full <- train.full[, -nzv]
        test.full <- test.full[, -nzv]
    }
    tooHigh <- findCorrelation(cor(train.full), .9)
    if (length(tooHigh) > 0) {
        train.full <- train.full[, -tooHigh]
        test.full <- test.full[, -tooHigh]
    }

    ## save the test set results in a data.frame
    if (length(test.outcome) > 0)
        testResults <- data.frame(obs=test.outcome)

    if (opt$debug) save.image('debug.Rdata')

    ## benchmark the specified models
    tuned.list <- list()
    accu <- data.frame()
    top.f <- data.frame()
    for (model in models) {
        if (opt$feature_selection) {
            ctrl <- rfeControl(method = "repeatedcv",
                               repeats = 5, number=10,
                               saveDetails = TRUE)
            ctrl$functions <- rfFuncs
            set.seed(721)
            tuned <- rfe(train.full,
                         train.outcome,
                         sizes = seq(10, ncol(train.full)-10, by=10),
                         metric = "RMSE",
                         ntree = 1000,
                         rfeControl = ctrl)


        } else {
            tuned <- regression.tune(train.full, train.outcome, model)
            ## if(is.na(tuned) | is.null(tuned)) next
            if (class(tuned) != 'train') {
                cat("Warning message:\nModel ", model, " failed.\n")
                next
            }
        }

        tuned.list[[model]] <- tuned

        if (opt$verbose) {
            print(tuned)
            print(accuracy(tuned))
        }

        ## add a new column - Model
        accu <- rbind(accu, cbind(tuned$resample, Model=tuned$method))

        if (length(test.outcome) > 0)
            testResults[model] <- predict(tuned, test.full)
        imp <- varImp(tuned)
        top.f <- rbind(top.f,
                       data.frame(imp$importance[order(imp$importance,
                                                       decreasing=T),,drop=FALSE],
                                  Model=model))
        ## plot top features
        pimp <- plot.imp(imp, tax.16s, main=paste(label, model))
        print(pimp, position=c(0, 0, 0.56, 1))

        save.image(sprintf("%s.Rdata", opt$output))
    }
    accu$Field <- label
    accuracies <- rbind(accuracies, accu)
    top.f$Field <- label
    top.features <- rbind(top.features, top.f)

    ## if (opt$verbose) print(accuracies)

    big.tuned.list[[label]] <- tuned.list

    if (length(tuned.list) > 1) {
        ## compare the model performances
        resamp <- resamples(tuned.list)
        m.diff <- diff(resamp)
        if (opt$verbose) print(summary(m.diff))
        print(dotplot(m.diff, main=label))
    }

    ## plot yhat vs. obs
    if (length(test.outcome) > 0) {
        method.names <- names(testResults)
        obs <- testResults[,1]
        for(i in 2:length(testResults)) {
            pred <- testResults[,i]
            plot(pred ~ obs, main = label,
                 xlab=method.names[1], ylab=method.names[i])
            abline(0, 1, col="red")
            ## mtext(paste(c("RMSE=", "R^2="),
            ##             c(RMSE())))
            rmse <- format(round(caret::RMSE(pred, obs), 2), nsmall=2)
            rsq <- format(round(caret::R2(pred, obs), 2), nsmall=2)
            legend("topleft", text.col="blue", "ab",
                   paste(c("RMSE","R^2 "), c(rmse, rsq), sep='=', collapse='\n'))
        }
    }
}

dev.off()

save.image(sprintf("%s.Rdata", opt$output))
