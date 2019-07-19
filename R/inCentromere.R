#' Return an index of segments that are in the centromere, for a given chromosome
#' 
#' @param start.pos start position of the segment
#' @param end.pos end position of the segment
#' @param ref.cent.start start position of the centromere
#' @param ref.cent.end end position of the centromere
#' 
#' @details This function returns a logical vector with value of TRUE if the corresponding segment crosses
#' the centromere.
#' 
#' @return logical index of segments crossing the centromere
#' 
#' @export

# ---------------------------------- in.centromere ------------------------------ #
# logical function to 
in.centromere <- function( start.pos, end.pos, ref.cent.start, ref.cent.end )
# input:
# start.pos (integer), start position of the segment (in bp)
# end.pos (integer), ending position of the segment
# ref.cent.start (integer), starting position of the centromere
# ref.cent.end (integer), ending position of the centromere
#
# output: index of logical values denoting segments in the centromere
#
{
  # remove any segments that cross the centromere
  return( start.pos < ref.cent.start & end.pos > ref.cent.end |
    start.pos > ref.cent.start & start.pos < ref.cent.end |
    end.pos > ref.cent.start & end.pos < ref.cent.end )
  
  
}
# ------------------------------------------------------------------------------- #
