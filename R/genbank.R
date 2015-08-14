
# Match a pattern in a fixed substring.
#
# posmatch(str, start, end, pattern) returns TRUE if pattern (a regexp) matches
# the substr(str, start, end).  Returns FALSE if pattern does not match or if
# str has fewer than end characters.
posmatch <- function(str, start, end, pattern) {
  if (nchar(str) < end) {
    return(FALSE)
  } else {
    return(str_detect(substr(str, start, end), pattern))
  }
}

# Returns TRUE if a character vector of lines appears to begin with a GenBank
# header.
has_header <- function(lines) {
  return(posmatch(lines[1], 21, 46, "Genetic Sequence Data Bank"))
}

# Parse the header section of a GenBank file.
#
# Input should be character vector of length 9; if longer, only the first 9
# lines are used.
#
# Returns a list containing the filename, date, release, title, loci, bases,
# and reports for the file.  If lines does not begin with a valid header, an 
# empty list is returned.
parse_header <- function(lines) {
  if (!has_header(lines)) {
    return(list())
  }
  
  return(list(
    # format says filename is 9 characters, but some appear longer?
    filename = str_trim(substr(lines[1], 1, 20)),
    date = str_trim(lines[2]),
    release = str_trim(substr(lines[4], 48, 52)),
    title = str_trim(lines[6]),
    loci = strtoi(substr(lines[8], 1, 8)),
    bases = strtoi(substr(lines[8], 16, 26)),
    reports = strtoi(substr(lines[8], 40, 47))
  ))
}

# Parse a group of lines representing a single keyword section.
parse_keyword_group <- function(group) {
  keys <- str_trim(substr(group, 1, 12))
  vals <- str_trim(substr(group, 13, nchar(group)))
  keyvals <- join_unkeyed_lines(keys, vals, sep="\n")
  keys <- keyvals$keys
  vals <- keyvals$values
  args <- vals[-1]
  names(args) <- keys[-1]
  args["value"] <- vals[1]
  return(list(key=keys[1], value=as.list(args)))
}

# Parse a origin character vector into a DNAString object.
parse_origin <- function(origin) {
  # use origin[-1] to remove keyword line
  origin[-1] %>%
    str_replace_all("\\s|\\d", "") %>%
    paste0(collapse="") %>%
    Biostrings::DNAString()
}

# Parse a feature character vector in a list of features.
parse_features <- function(features) {
  groups <- split_by_indent(features[-1])  # -1 to remove the keyword line
  
  parse_feature <- function(lines) {
    header <- lines[1]
    lines <- substr(lines[-1], 22, nchar(lines[-1]))  # remove leading ws
    key <- str_trim(substr(header, 6, 21))
    location_str <- str_trim(substr(header, 22, nchar(header)))
    location <- parse_location(location_str)
    keys <- str_match(lines, "^/(.+)=")[ ,2]
    vals <- str_match(lines, "^/.+=(.+)")[ ,2]
    vals[is.na(vals)] <- lines[is.na(vals)]  # use whole line when continuation
    keyvals <- join_unkeyed_lines(keys, vals, sep="")
    ret <- keyvals$values %>% 
      str_replace("^\"", "") %>% 
      str_replace("\"$", "") %>%
      as.list()
    names(ret) <- keyvals$keys
    ret$key <- key
    ret$location <- location
    ret$location_str <- location_str
    return(list(key=key, value=ret))
  }
  
  return(group_lists_by_key(lapply(groups, parse_feature)))
}

parse_genbank <- function(file) {
  lines <- readLines(file)
  
  if (has_header(lines)) {
    header <- parse_header(lines[1:9])  # header is only first 9 lines
    lines <- lines[-(1:9)]  # strip off the header lines
  } else {
    header <- NULL
  }
  
  groups <- split_by_indent(lines)
  names(groups) <- lapply(groups, f(x, extract_keywords(x, first_only=T)))
  features <- groups$FEATURES
  groups$FEATURES <- NULL
  origin <- groups$ORIGIN
  groups$ORIGIN <- NULL
  if ("//" %in% names(groups)) {
    # remove end of file character
    groups[["//"]] <- NULL
  }
  
  names(groups) <- NULL
  gbk <- group_lists_by_key(lapply(groups, parse_keyword_group))
  # singlets appear only once per file and have no subkeys
  singlets <- c("KEYWORDS", "VERSION", "ACCESSION", "DEFINITION", "LOCUS")
  for (singlet in singlets) {
    gbk[[singlet]] <- gbk[[singlet]][[1]]$value
  }
  # onesies appear only once per file but have subkeys
  onesies <- c("SOURCE")
  for (onesie in onesies) {
    gbk[[onesie]] <- gbk[[onesie]][[1]]
  }
  
  gbk$header <- header
  gbk$sequence <- parse_origin(origin)
  gbk$features <- parse_features(features)
  
  return(gbk)
}
