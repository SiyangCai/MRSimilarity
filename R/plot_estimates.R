#' Construct a scatter plot of instrument effects and causal effect estimates.
#'
#' Plot instrument-exposure effects and instrument-outcome effects (with corresponding 95% confidence intervals). Overlay causal effect estimates as slopes to allow visual inspection of similarity between causal effect estimates.
#' @param beta.exposure A vector of instrument-exposure association estimates.
#' @param beta.outcome A vector of instrument-outcome association estimates from a common instrument set.
#' @param se.exposure A vector of standard error estimates for the set of instrument-exposure association estimates.
#' @param se.outcome A vector of stand error estimates for the set of instrument-outcome association estimates.
#' @param b A vector of slope estimates for a given exposure-outcome pair.
#' @param beta.intercepts A vector of intercepts corresponding to the set of causal effect estimates given. Only the MR-Egger intercept should be non-zero.
#' @param title_text A string object specifying the title of the customised coverage plot.
#' @param xlab A string object specifying the x-axis title for the scatter plot.
#' @param ylab A string object specifying the y-axis title for the scatter plot.
#' @param methods A vector detailing the name of chosen causal effect estimators.
#' @param xsize A numeric object specifying the desired size of the x-axis tick markers. The default value is 22.
#' @param ysize A numeric object specifying the desired size of the y-axis tick markers. The default value is 22.
#' @param titlesize A numeric object specifying the desired size of the plot title. The default value is 24.
#' @param opacity A numeric object specifying the opacity of the scatter plot markers (between 0 and 1). The default value is 1.
#'
#' @return Returns a plotly plot object.
#' @export
#'
#' @examples SNP.plot <- MR_2s_scatter_auto(beta.exposure = SNP_LDL.beta, beta.outcome = SNP_CHD.beta, se.exposure = SNP_LDL.beta.se, se.outcome = SNP_CHD.beta.se, b = ests.beta, beta.intercepts = c(0, -0.01, 0, 0), title.text = "Causal effect estimates of LDL cholesterol on coronary heart disease risk", xlab = "SNP-LDL effect", ylab = "SNP-CHD effect", methods = c("IVW", "MRE", "WM", "MBE"), xsize = 25, ysize = 25, titlesize = 25)

MR_scatter_estimates <-
  function(beta.exposure,
           beta.outcome,
           se.exposure,
           se.outcome,
           b,
           beta.intercepts,
           title_text,
           xlab,
           ylab,
           methods,
           xsize = 22,
           ysize = 22,
           titlesize = 24,
           opacity = 1) {

    # color vector - enough to match with each 18 TwoSampleMR method if requested
    colors.vec <-
      c(
        "coral",
        "dodgerblue",
        "darkgreen",
        "blueviolet",
        "purple",
        "brown",
        "pink",
        "DeepSkyBlue",
        "Aquamarine",
        "Maroon",
        "LightSlateGray",
        "Aqua",
        "YellowGreen",
        "Magenta",
        "LemonChiffon",
        "DarkKhaki",
        "Thistle",
        "DarkSlateBlue"
      )


    # set up plot object with individual instrument effects
    p <- plotly::plot_ly(type = "scatter", mode = "markers") %>%
      plotly::add_trace(
        x = beta.exposure,
        y = beta.outcome,
        marker = list(color = "black", opacity = opacity),
        name = "Genetic Instruments",
        error_x = list(
          array = se.exposure,
          color = 'rgba(192, 192, 192, 0.3)'
          #color = "gray",
          #opacity = 0.15
        ),
        error_y = list(
          array = se.outcome,
          color = 'rgba(192, 192, 192, 0.3)'
          #color = "gray",
          #opacity = opacity
        )
      ) %>%
      plotly::layout(
        legend = list(
          x = 1.02,
          y = 0.50,
          font = list(size = 16),
          bordercolor = "gray",
          borderwidth = 2
        ),
        title = list(
          text = title_text,
          x = 0.53,
          y = 0.94,
          anchor = "center",
          yanchor = "top",
          font = list(size = titlesize)
        ),
        xaxis = list(
          tickfont = list(size = 17),
          title = list(
            text = xlab,
            font = list(size = xsize),
            standoff = 10
          ),
          y = 0.1,
          standoff = 50
        ),
        yaxis = list(
          tickfont = list(size = 17),
          title = list(text = ylab, font = list(size = ysize))
        ),
        margin = list(r = 25, t = 50)
      )

    # iterate through effect data and add each slope
    for (i in c(1:length(methods))) {
      p <- p %>% plotly::add_trace(
        x = c(0, max(beta.exposure) * 1.1),
        y = c(beta.intercepts[i],
              max(beta.exposure) * b[i] * 1.1),
        name = methods[i],
        mode = "lines",
        line = list(
          width = 2,
          dash = "solid",
          color = colors.vec[i]
        )
      )
    }

    return(p)
  }
