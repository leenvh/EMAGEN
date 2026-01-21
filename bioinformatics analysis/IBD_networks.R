## script adapted from https://github.com/vivaxgen

library(igraph)
library(unikn)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

mytriangle <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1 / 200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }

  symbols(
    x = coords[, 1], y = coords[, 2], bg = vertex.color,
    stars = cbind(vertex.size, vertex.size, vertex.size),
    add = TRUE, inches = FALSE
  )
}

add_shape("triangle",
  clip = shapes("circle")$clip,
  plot = mytriangle
)


mystar <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1 / 200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  norays <- params("vertex", "norays")
  if (length(norays) != 1 && !is.null(v)) {
    norays <- norays[v]
  }

  mapply(coords[, 1], coords[, 2], vertex.color, vertex.size, norays,
    FUN = function(x, y, bg, size, nor) {
      symbols(
        x = x, y = y, bg = bg,
        stars = matrix(c(size, size / 2), nrow = 1, ncol = nor * 2),
        add = TRUE, inches = FALSE
      )
    }
  )
}

add_shape("star",
  clip = shape_noclip,
  plot = mystar, parameters = list(vertex.norays = 5)
)


.reorder.metadata <-
  function(IBD, metadata, sample.col = "Sample") {
    samples <- unique(c(IBD[, "sample1"], IBD[, "sample2"]))
    sample.order <- match(metadata[, sample.col], samples)
    actual.length <- length(na.omit(sample.order))
    
    metadata <- metadata[order(sample.order), ]
    metadata[1:actual.length, ]
  }


.get.labels <-
  function(metadata,
           label.col,
           label.palette = NULL) {
    labels <- metadata[, label.col]
    labels[is.na(labels)] <- NaN
    
    labels <- factor(labels)
    
    if (is.null(label.palette)) {
      labels
    } else {
      factor(labels, levels = names(label.palette))
    }
  }


.generate.label.palette <- function(labels, palette = "custom") {
  label.names <- levels(labels)
  num.labels <- length(label.names)
  
  if (palette == "custom") {
    custom_palette <- c(
      "#A6CEE3", 
      "#238b8d", 
      "#a72e61", 
      "#FDBF6F", 
      "#CAB2D6", 
      "#FFFF99", 
      "#1F78B4", 
      "#33A02C", 
      "#E31A1C", 
      "#FF7F00", 
      "#6A3D9A", 
      "#B15928",
      "#F4A582", 
      "#92C5DE", 
      "#D9F0D3"  
    )
    
    if (num.labels > length(custom_palette)) {
      label.palette <- rep(custom_palette, length.out = num.labels)
    } else {
      label.palette <- custom_palette[1:num.labels]
    }
  } else if (palette == "unikn") {
    label.palette <- unikn::usecol(pal = pal_unikn, n = num.labels)
  } else {
    label.palette <- palette.colors(num.labels, palette = palette)
  }
  
  names(label.palette) <- label.names
  label.palette
}



.generate.label.colours <-
  function(labels, label.palette)
    label.palette[labels]


.create.edgelist <-
  function(IBD)
    IBD[, c("sample1", "sample2", "fract_sites_IBD")]


.create.vertices <-
  function(IBD,
           metadata = NULL,
           sample.col = "Sample") {
    if (is.null(metadata)) {
      samples <- unique(c(IBD[, "sample1"], IBD[, "sample2"]))
      vertices <- data.frame(samples)
      
    } else {
      metadata <-
        .reorder.metadata(IBD, metadata, sample.col = sample.col)
      vertices <- data.frame(metadata[, sample.col])
    }
    
    names(vertices) <- sample.col
    
    vertices
  }


.plot.IBD <-
  function(IBD,
           IBD.cutoff,
           metadata,
           label.col,
           label.palette = NULL,
           coords.cache = NULL,
           unlabelled = TRUE) {
    edgelist <- .create.edgelist(IBD)
    metadata <- .reorder.metadata(IBD, metadata)
    vertices <- .create.vertices(IBD, metadata = metadata)
    
    
    d <- edgelist[edgelist[, "fract_sites_IBD"] >= IBD.cutoff, ]
    IBD.graph <-
      graph_from_data_frame(d, directed = FALSE, vertices = vertices)
    
    
    # store metadata about the IBD graph
    if (!is.null(coords.cache)) {
      coords <- coords.cache[[as.character(IBD.cutoff)]]
    } else {
      coords <- layout_(IBD.graph, nicely())
    }
    
    if (is.null(coords)) {
      coords <- layout_(IBD.graph, nicely())
      coords.cache[[as.character(IBD.cutoff)]] <- coords
    }
    
    
    # add label and colour to the nodes
    if (is.null(label.palette)) {
      labels <- .get.labels(metadata, label.col)
      label.palette <- .generate.label.palette(labels)
      
    } else {
      labels <-
        .get.labels(metadata, label.col, label.palette = label.palette)
    }
    
    label.colours <-
      .generate.label.colours(labels, label.palette)
    
    IBD.graph <-
      set_vertex_attr(IBD.graph, "color", value = label.colours)
    
    # plot the IBD graph
    percent.cutoff <- round(IBD.cutoff * 100, digits = 2)
    
    if (unlabelled) {
      plot(
        IBD.graph,
        layout = coords,
        vertex.size = 6,
        vertex.label = NA,
        main = paste0("IBD >=", percent.cutoff, "%")
      )
    } else {
      plot(
        IBD.graph,
        layout = coords,
        vertex.size = 6,
        vertex.label.cex = 0.3,
        main = paste0("IBD >=", percent.cutoff, "%")
      )
    }
    
    
    coords.cache
  }

plot.IBD <-
  function(file,
           IBD,
           IBD.cutoffs,
           metadata,
           label.col,
           label.palette = NULL,
           onefile = TRUE,
           coords.cache = NULL,
           unlabelled = TRUE) {
    # extract relevant data to create a graph of IBDs
    edgelist <- .create.edgelist(IBD)
    metadata <- .reorder.metadata(IBD, metadata)
    
    # add label and colour to the nodes
    if (is.null(label.palette)) {
      labels <- .get.labels(metadata, label.col)
      label.palette <- .generate.label.palette(labels)
      
    } else {
      labels <-
        .get.labels(metadata, label.col, label.palette = label.palette)
    }
    
    
    if (onefile) {
      sqrt.n.plots <- ceiling(sqrt(length(IBD.cutoffs)))
      bottom.left.index <- (sqrt.n.plots - 1) * sqrt.n.plots + 1
      size <- 6 * sqrt.n.plots
      
      pdf(
        file = file,
        width = size,
        height = size
      )
      
      par(mar = c(1, 11, 1, 1) + 0.1,
          mfrow = c(sqrt.n.plots, sqrt.n.plots),
          cex.main = 1.8)
    }
    
    for (i in seq_along(IBD.cutoffs)) {
      IBD.cutoff <- IBD.cutoffs[i]
      
      if (!onefile) {
        newfile <- paste0(file, "_IBD", IBD.cutoff, ".pdf")
        pdf(file = newfile)
        par(mar = c(1, 6, 1, 1) + 0.1)
      }
      
      coords.cache <- .plot.IBD(
        IBD,
        IBD.cutoff,
        metadata,
        label.col,
        label.palette = label.palette,
        coords.cache = coords.cache,
        unlabelled = unlabelled
      )
      
      if (!onefile) {
        # add the legend about the IBD graph
        legend(
          "bottomleft",
          legend = levels(labels),
          fill = label.palette,
          cex = 0.8,
          title = label.col,
          inset = c(-0.2, 0),
          xpd = NA
        )
        dev.off()
      } else {
        if (i == bottom.left.index) {
          # add the legend about the IBD graph
          legend(
            "bottomleft",
            legend = levels(labels),
            fill = label.palette,
            cex = 1.5,
            title = label.col,
            inset = c(-0.3, 0),
            xpd = NA
          )
        }
      }
    }
    
    if (onefile)
      dev.off()
    
    invisible(coords.cache)
  }


# include -5% IBD
IBD.cutoffs <- c(0.00390625, 0.0078125, 0.015625, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1) * 0.95

metadata <-
  read.csv("../District_and_Sample_ID_Metadata_oct11.csv")


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


metadata[, "district"] <- metadata[, "District"]

label.col <- "District"
IBD.plot.file <- "IBD_district.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)



######## R622I ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the Pfk13_R622I column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfk13_R622I")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)



metadata[, "K13_R622I"] <- metadata[,"Pfk13_R622I"]

label.col <- "Pfk13_R622I"
IBD.plot.file <- "IBD_Pfk13_R622I.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)




######## A675V ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the A675V column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfk13_A675V")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)



metadata[, "K13_A675V"] <- metadata[,"Pfk13_A675V"]

label.col <- "Pfk13_A675V"
IBD.plot.file <- "IBD_Pfk13_A675V.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)


######## P441L ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the P441L column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfk13_P441L")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)



metadata[, "K13_P441L"] <- metadata[,"Pfk13_P441L"]

label.col <- "Pfk13_P441L"
IBD.plot.file <- "IBD_Pfk13_P441L.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)

######## CVIET ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")
metadata_mutations$CVIET<-"No"
metadata_mutations$CVIET[metadata_mutations$Pfcrt_K76T=="Yes" & metadata_mutations$Pfcrt_M74I=="Yes" & metadata_mutations$Pfcrt_N75E=="Yes"] <- "Yes"


# Add the CVIET column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "CVIET")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


label.col <- "CVIET"
IBD.plot.file <- "IBD_CVIET.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)



######## K76T ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the K76T column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfcrt_K76T")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


label.col <- "Pfcrt_K76T"
IBD.plot.file <- "IBD_Pfcrt_K76T.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)


######## Pfcrt_A220S ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the Pfcrt_A220S column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfcrt_A220S")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


label.col <- "Pfcrt_A220S"
IBD.plot.file <- "IBD_Pfcrt_A220S.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)


######## Pfcrt_I356T ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the Pfcrt_I356T column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfcrt_I356T")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


label.col <- "Pfcrt_I356T"
IBD.plot.file <- "IBD_Pfcrt_I356T.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)



######## Pfcrt_M74I ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the Pfcrt_M74I column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfcrt_M74I")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


label.col <- "Pfcrt_M74I"
IBD.plot.file <- "IBD_Pfcrt_M74I.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)



######## Pfcrt_N75E ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the Pfcrt_N75E column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfcrt_N75E")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


label.col <- "Pfcrt_N75E"
IBD.plot.file <- "IBD_Pfcrt_N75E.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)


######## Pfcrt_Q271E ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the Pfcrt_N75E column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfcrt_Q271E")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


label.col <- "Pfcrt_Q271E"
IBD.plot.file <- "IBD_Pfcrt_Q271E.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)

######## Pfcrt_N326S ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the Pfcrt_N326S column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfcrt_N326S")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


label.col <- "Pfcrt_N326S"
IBD.plot.file <- "IBD_Pfcrt_N326S.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)

######## Pfcrt_R371I ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the Pfcrt_R371I column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfcrt_R371I")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


label.col <- "Pfcrt_R371I"
IBD.plot.file <- "IBD_Pfcrt_R371I.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)


######## Pfdhfr_C59R ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the Pfdhfr_C59R column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfdhfr_C59R")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


label.col <- "Pfdhfr_C59R"
IBD.plot.file <- "IBD_Pfdhfr_C59R.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)



######## Pfdhfr_N51I ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the Pfdhfr_N51I column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfdhfr_N51I")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


label.col <- "Pfdhfr_N51I"
IBD.plot.file <- "IBD_Pfdhfr_N51I.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)

######## Pfdhfr_S108N ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the Pfdhfr_S108N column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfdhfr_S108N")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


label.col <- "Pfdhfr_S108N"
IBD.plot.file <- "IBD_Pfdhfr_S108N.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)



######## Pfdhps_A437G ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the Pfdhps_A437G column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfdhps_A437G")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


label.col <- "Pfdhps_A437G"
IBD.plot.file <- "IBD_Pfdhps_A437G.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)


######## Pfdhps_S436A ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the Pfdhps_S436A column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfdhps_S436A")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


label.col <- "Pfdhps_S436A"
IBD.plot.file <- "IBD_Pfdhps_S436A.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)


######## Pfdhps_K540E ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the Pfdhps_K540E column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfdhps_K540E")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


label.col <- "Pfdhps_K540E"
IBD.plot.file <- "IBD_Pfdhps_K540E.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)


######## Pfdhps_A581G ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the Pfdhps_A581G column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfdhps_A581G")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


label.col <- "Pfdhps_A581G"
IBD.plot.file <- "IBD_Pfdhps_A581G.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)


######## Pfmdr1_N86Y ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the Pfmdr1_N86Y column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfmdr1_N86Y")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


label.col <- "Pfmdr1_N86Y"
IBD.plot.file <- "IBD_Pfmdr1_N86Y.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)


######## Pfmdr1_Y184F ########

metadata_mutations <-
  read.csv("../Metadata_mutations.csv")

# Add the Pfmdr1_Y184F column to metadata
metadata <- merge(
  metadata, 
  metadata_mutations[, c("Sample_ID", "Pfmdr1_Y184F")], 
  by = "Sample_ID", 
  all.x = TRUE
)


# inefficient, but easiest
base <- metadata
for (i in 1:4) {
  copy <- base
  copy[, "Sample_ID"] <- paste0(copy[, "Sample_ID"], "-", i)
  metadata <- rbind(metadata, copy)
}

metadata[, "Sample"] <- metadata[, "Sample_ID"]

coords.cache <- vector(mode = "list", length = length(IBD.cutoffs))
names(coords.cache) <- as.character(IBD.cutoffs)

IBD.file <- "hmmIBD_Ethiopia_maf0.01_out.hmm_fract.txt"
IBD <- read.delim(IBD.file)


label.col <- "Pfmdr1_Y184F"
IBD.plot.file <- "IBD_Pfmdr1_Y184F.pdf"

coords.cache <- plot.IBD(
  IBD.plot.file,
  IBD,
  IBD.cutoffs,
  metadata,
  label.col,
  coords.cache = coords.cache
)



#metadata[grep("_D0$", metadata[, "external_id"], perl = TRUE), "Status"] <- "Day 0"
#metadata[grep("_D0?RP.$", metadata[, "external_id"], perl = TRUE), "Status"] <- "Recurrence"

#label.col <- "Status"
#IBD.plot.file <- "/home/ashley/Documents/Ethiopia/pv_IMPROV/DEploid/BEST/plots/filtered_non-swga_AF_monoclonals.DEploid.major_IBD_status.pdf"

#plot.IBD(IBD.plot.file,
#         IBD,
#         IBD.cutoffs,
#         metadata,
#         label.col,
#         coords.cache = coords.cache)