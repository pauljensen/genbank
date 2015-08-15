

pattern_idx <- "(\\d+)"
pattern_loc <- paste0("([<>])?", pattern_idx)
fullstr <- function(...) paste0("^", ..., "$")
pattern_single <- fullstr(pattern_idx)
pattern_site <- fullstr(pattern_loc, "\\^", pattern_loc)
pattern_span <- fullstr(pattern_loc, "\\.\\.", pattern_loc)
pattern_range <- fullstr(pattern_loc, "\\.", pattern_loc)
pattern_op <- fullstr("(\\w+)\\((.+)\\)")

simple_op <- function(operator, start="", end="", start_open="", end_open="") {
  list(operator=operator, 
       start=as.integer(start), 
       end=as.integer(end), 
       start_open=ifelse(nchar(start_open) > 0, start_open, NA), 
       end_open=ifelse(nchar(end_open) > 0, end_open,  NA))
}

parse_location <- function(l) {
  if (str_detect(l, pattern_single)) {
    # single base
    m <- str_match(l, pattern_single)
    return(simple_op("base", start=m[1,2], end=m[1,2]))
    
  } else if (str_detect(l, pattern_site)) {
    # site between bases
    m <- str_match(l, pattern_site)
    return(simple_op("site", start_open=m[1,2], start=m[1,3], 
                     end_open=m[1,4], end=m[1,5]))
    
  } else if (str_detect(l, pattern_span)) {
    # base span
    m <- str_match(l, pattern_span)
    return(simple_op("span", start_open=m[1,2], start=m[1,3], 
                     end_open=m[1,4], end=m[1,5]))
    
  } else if (str_detect(l, pattern_range)) {
    # base range
    m <- str_match(l, pattern_range)
    return(simple_op("range", start_open=m[1,2], start=m[1,3], 
                     end_open=m[1,4], end=m[1,5]))
    
  } else if (str_detect(l, pattern_op)) {
    # operator
    m <- str_match(l, pattern_op)
    arg_strs <- str_split(m[1,3], ",\\s*")[[1]]
    return(list(operator=m[1,2], args=lapply(arg_strs, parse_location)))
  }
}

extract_sequence <- function(loc, dna) {
  if (is.character(loc)) {
    loc <- parse_location(loc)
  }
  extract <- f(x, extract_sequence(x, dna))  # for recursion
  
  if (loc$operator == "base") {
    subseq(dna, loc$start, loc$end)
  } else if (loc$operator == "site") {
    DNAString("")
  } else if (loc$operator == "span") {
    subseq(dna, loc$start, loc$end)
  } else if (loc$operator == "complement") {
    Biostrings::complement(extract(loc$args[[1]])) # complement takes single arg
  } else if (loc$operator == "join") {
    do.call(c, lapply(loc$args, extract))
  } else if (loc$operator == "order") {
    lapply(loc$args, extract)
  } else {
    DNAString("")
  }
}
