#' Data from Experiment 1 (Wagenmakers, Ratcliff, Gomez, & McKoon, 2008)
#' 
#' Responses and response times from an experiment in which instruction manipulated speed and accuracy between blocks. This data was also analyzed by Heathcote and Love (2012) who were the first to use the 17 participants also included here. 
#' 
#' @details The data excludes the practie blocks but includes all trials. Variable \code{censor} can be used for excluding all trials also excluded from the papers using it namely uninterpretable response, too fast response (<180 ms), too slow response (>3 sec). Heathcote and Love (2012, p. 7) describe the data as follows:
#' 
#' We fit the LBA and LNR models to data from Wagenmaker et al.'s (2008) experiment one, where participants made decisions about whether a string of letters constituted a word. These lexical decisions were made about four types of stimuli, non-words (nw) and high-frequency (hf), low-frequency (lf), and very low-frequency (vlf) words. Participants made decisions either under speed or accuracy emphasis instructions in different experimental blocks. Accuracy blocks were preceded by the message "Try to respond accurately" and "ERROR" was displayed after each wrong response. Speed blocks were preceded by the message "Try to respond accurately" and "TOO SLOW" was displayed after each response slower than 0.75 s.We report analyses of data from 17 participants (31,412 data points) in their Experiment 1, including the 15 participants analyzed in Wagenmakers et al. (2008) and two extras (we thank Eric-Jan Wagenmakers for supplying this data).
#' 
#' @docType data
#' @keywords dataset
#' @name speed_acc
#' @usage speed_acc
#' @format A \code{data.frame} with 31,522 obs. and 9 variables:
#' \describe{
#'  \item{id}{participant id}
#'  \item{block}{block number}
#'  \item{condition}{\code{accuracy} for blocks with accuracy instructions; \code{speed} for blocks with speed instruction}
#'  \item{stim}{unique identifier of stimulus, stimuli are nested in frequency conditions}
#'  \item{stim_cat}{category of stimulus, either word or non-word}
#'  \item{frequency}{"high frequency word", "low frequency word", "very low frequency word", or non-words derived from the first three categories}
#'  \item{response}{\code{word}, \code{nonword}, or not interpretable response (\code{error}, i.e., pushed a button, but not the right one and also not the one next to the right button)}
#'  \item{rt}{response time in seconds}
#'  \item{censor}{boolean indicating whether or not a response should be eliminated prior to analysis; ninterpretable response, too fast response (<180 ms), too slow response (>3 sec)}
#' }
#'
#' @source Wagenmakers, E.-J., Ratcliff, R., Gomez, P., & McKoon, G. (2008). A diffusion model account of criterion shifts in the lexical decision task. \emph{Journal of Memory and Language}, 58(1), 140-159.
#' 
#' @references Heathcote, A., & Love, J. (2012). Linear deterministic accumulator models of simple choice. \emph{Frontiers in Psychology: Cognitive Science}, 3, 292. doi:10.3389/fpsyg.2012.00292
#' 
#' @example examples/examples.speed_acc.R
#'
NULL