
f <- pryr::f

# Split lines into a set of groups where items 2+ are indented more than the
# first item.
#
# Input is a character vector.  Returns a list of character vectors.
# Inputting a character vector of length 0 returns an empty list.
split_by_indent <- function(lines) {
  if (length(lines) == 0) {
    return(list())
  } else if (length(lines) == 1) {
    return(list(lines))
  }
  
  depth <- nchar(str_match(lines, "^\\s*")[ ,1])
  groups <- list()
  group_idx <- 0
  current_start <- 1
  for (i in 2:length(lines)) {
    if (depth[i] == depth[current_start]) {
      # starting new group
      group_idx <- group_idx + 1
      groups[[group_idx]] <- lines[current_start:(i-1)]
      current_start <- i
    }
  }
  # push the last group
  groups[[group_idx+1]] <- lines[current_start:length(lines)]
  return(groups)
}

# Returns GenBank keywords and subkeywords in columns 1-12 given a character
# vector.  If first_only, return only the keyword from the first line.
extract_keywords <- function(lines, first_only=F) {
  keywords <- str_trim(substr(lines, 1, 12))
  if (first_only) {
    return(keywords[1])
  } else {
    return(keywords)
  }
}

# Concatenates values without a key.  Inputs are two character vectors (keys
# and values), where values with key == "" or NA are concatenated onto the 
# previous value with separator sep.  Returns a list of character vectors with 
# names "keys" and "values".
#
# Example:
#   keys   = ["a",      "",       "b"]
#   values = ["Line 1", "Line 2", "Line 3"]
# becomes
#   list(keys=["a", "b"], values=["Line 1\nLine2", "Line 3"])
join_unkeyed_lines <- function(keys, values, sep="\n") {
  keys[is.na(keys)] <- ""
  nkeys <- sum(sapply(keys, f(nchar(x) > 0)))
  new_keys <- character(nkeys)
  new_values <- character(nkeys)
  current_key <- 0
  for (i in seq_along(keys)) {
    if (nchar(keys[i]) > 0) {
      current_key <- current_key + 1
      new_keys[current_key] <- keys[i]
      new_values[current_key] <- values[i]
    } else {
      new_values[current_key] <- paste(new_values[current_key], values[i], 
                                       sep=sep)
    }
  }
  return(list(keys=new_keys, values=new_values))
}

# Combine objects in a list based on keys.
# Given a list of lists (each with names "key" and "value"), consolidate items
# with identical keys.  lists should not have names.
#
# keyfun and valfun can be used for custom key and value generation.
#
# Example:
#   list( list(key=k1, value=v1),
#         list(key=k2, value=v2),
#         list(key=k1, value=v3)  )
# becomes
#   list( k1=list(v1, v3),
#         k2=list(v2)      )
group_lists_by_key <- function(lists, keyfun=f(x$key), valfun=f(x$value)) {
  keys <- sapply(lists, keyfun)
  vals <- lapply(lists, valfun)
  grouped <- list()
  for (key in unique(keys)) {
    grouped[[key]] <- vals[key == keys]
  }
  return(grouped)
}

