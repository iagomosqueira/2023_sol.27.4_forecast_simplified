# utilities.R - DESC
# /utilities.R

# Copyright Iago MOSQUEIRA (WMR), 2021
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(ftExtra)

# interimTable {{{

interimTable <- function(run, fy=dims(run)$maxyear) {

  ay <- ac(fy - 2)
  iy <- ac(fy - 1)

  tabs <- FLQuants(
  `F~ages\\ 2-6~ (ay)`=fbar(run)[, ay],
  `SSB (iy)`=ssb(run)[, iy],
  `R~age\\ 1~ (ay, iy)`=rec(run)[, ay],
  `Total catch (ay)`=catch(run)[, ay],
  `Projected landings (ay)`=landings(run)[, ay],
  `Projected discards (ay)`=discards(run)[, ay])

  names(tabs) <- gsub("iy", iy, names(tabs))
  names(tabs) <- gsub("ay", ay, names(tabs))

  tab <- lapply(tabs, c)
  tab <- data.frame(Variable=names(tab), Value=unlist(tab), row.names=NULL)

  # icesRound 1st row
  tab[1, "Value"] <- icesRound(as.numeric(tab[1, "Value"]))

  # ROUND others
  tab[2:6, "Value"] <- round(as.numeric(tab[2:6, "Value"]))

  return(tab)

} # }}}

# catchOptionsTable {{{

catchOptionsTable <- function(runs, advice, tac=advice, ages, discards.ages) {

  sy <- dims(runs[[1]])$maxyear
  fy <- ac(sy - 1)
  iy <- ac(sy - 2)
  sy <- ac(sy)

  # SET catch options table rows
  rows <- c(
    advice = "MSY approach",
    Fmsy = "F~MSY~",
    lFmsy = "F = F~MSY lower~",
    lFmsyadvice = "F = F~MSY lower~ advice",
    uFmsy = "F = F~MSY upper~",
    F0 = "F = 0",
    Fpa = "F~pa~",
    F05noAR = "F~p.05~ without AR",
    Flim = "F~lim~",
    Bpa = paste0("SSB (", sy, ") = B~pa~"),
    Blim = paste0("SSB(", sy, ")=B~lim~"),
    MSYBtrigger = paste0("SSB(", sy, ") = MSY B~trigger~"),
    lastF = paste0("F~", iy, "~"),
    rotac = "Roll-over TAC")

  tabs <- lapply(runs[names(rows)], function(x) {
    model.frame(FLQuants(
      # 
      c(metrics(window(x, start=fy, end=fy),
        list(catch=catch, wanted=landings, unwanted=discards, F=fbar,
        Fwanted=function(z) Fwanted(z, ages),
        Funwanted=function(z) Funwanted(z, discards.ages))),
      metrics(window(x, start=sy, end=sy),
        list(SSB=ssb)))), drop=TRUE)
    }
  )

  # COMBINE as single table

  tab <- rbindlist(tabs, idcol="basis")

  # ADD SSB change, 100 - (old / new) * 100
  ssbfy <- c(ssb(runs[[1]])[, ac(fy)])
  tab[, ssbchange:=(SSB - ssbfy) / ssbfy * 100]

  # ADD TAC change
  tab[, tacchange:=(catch - c(tac)) / c(tac) * 100]

  # ADD advice change
  tab[, advicechange:=(catch - c(advice)) / c(advice) * 100]

  # FIX zeroes
  tab[, (2:11) := lapply(.SD, function(x)
    ifelse(abs(x) < 0.0001, 0, x)), .SDcols=2:11]

  tab$basis <- rows

  # CALL round and icesRound
  tab[, c('catch', 'wanted', 'unwanted', 'SSB') := lapply(.SD, round,
    digits=0), .SDcols = c('catch', 'wanted', 'unwanted', 'SSB')]
  tab[, c('ssbchange', 'tacchange', 'advicechange') := lapply(.SD, round,
    digits=2), .SDcols = c('ssbchange', 'tacchange', 'advicechange')]
  tab[, c('F', 'Fwanted', 'Funwanted') := lapply(.SD, icesRound),
      .SDcols = c('F', 'Fwanted', 'Funwanted')]

  cnms <- c("Basis", "Total catch (FY)", "Projected landings (FY)",
    "Projected discards (FY)", "F~total~ (ages AGES) (FY)",
    "F~projected\\ landings~ (ages AGES) (FY)",
    "F~projected\\ discards~ (ages DISC) (FY)", "SSB (SY)", "% SSB change",
    "% TAC change", "% advice change")

    cnms <- gsub("FY", fy, cnms)
    cnms <- gsub("SY", sy, cnms)
    cnms <- gsub("AGES", paste0(ages[1], "-", ages[2]), cnms)
    cnms <- gsub("DISC", paste0(discards.ages[1], "-", discards.ages[2]), cnms)

  setnames(tab, cnms)

  return(tab)
}

# }}}

# flextableCatchOptions {{{

flextableCatchOptions  <- function(x) {

  # ADD rows for sections
  x <- rbind(x[NA], x[1,], x[NA], x[-1])

  x[1, 1] <- "ICES advice basis"
  x[3, 1] <- "Other scenarios"

  # CONVERT to flextable
  ft <- flextable(x)

  # ADD footnotes
  ft <- footnote( ft, i = 1, j = c(2, 4, 5, 9, 11),
  value = as_paragraph(c(
    "Differences between the total catch and the sum of projected landings and discards result from rounding.",
    "Including BMS landings, assuming average estimated discard rate by age 2020-2022.",
    "SSB 2025 relative to SSB 2024.",
    "Total catch in 2024 relative to the advice value 2023 and TAC (both 9152 tonnes).",
    "Bpa and MSY Btrigger cannot be achieved in 2025, even with zero catches."
    )),
  ref_symbols = c("*", "**", "^", "^^", "##"),
  part = "header")

  # ADD footnote with as_sub
  ft <- footnote(ft, i = 1, j = 10,
  value = as_paragraph("F", as_sub("projected landings"), " and F", as_sub("projected discards"), " do not add up to F", as_sub("total"), " as they are calculated using different ages."),
  ref_symbols=c("#"), part="header")

  # MERGE section rows
  ft <- merge_h(ft, i = 1, part = "body")
  ft <- merge_h(ft, i = 3, part = "body")

  # SET vertical alignment
  ft <- valign(ft, valign = "bottom", part = "header")

  # SET markdown columns
  ft <- colformat_md(ft, j="Basis")
  ft <- colformat_md(ft, j=5:7, part="header")

  # ft <- autofit(ft, add_w = 4, part = c("body"),
  #   unit = "cm", hspans = "none")
  ft <- width(ft, 1:11, (c(3, rep(2, 10)) / 23) * 21, unit="cm")

  return(ft)
}
# }}}

# flextable(FLPar) {{{
setAs('FLPar', 'flextable',
  function(from) {

    # LIST
    dat <- as.list(from)
    names(dat) <- dimnames(from)$params

    # ROUND
    dat[c(2,5,6,7,8,9,10)] <- lapply(dat[c(2,5,6,7,8,9,10)], icesRound)
    dat[c(1,3,4)] <- lapply(dat[c(1,3,4)], round)

    #
    ft <- autofit(flextable(as.data.frame(dat)))

    return(ft)
  }
)
# }}}

# --- COMPARE runs {{{

compareFLStocks <- function(new, old, what="stock.n") {

  fy <- dims(old)$maxyear
  yrs <- ac(seq(fy - dims(old)$age + 1, fy))
  ags <- dimnames(old)$age

  dif <- do.call(what, list(new))[, yrs] / do.call(what, list(old))[, yrs]

  tab <- cbind(Age=ags, as.data.frame(format(round(dif[drop=TRUE], 2),
    digits=2)))

  return(tab)
}
# }}}
