# model.R - Run forecast, using Mohn's rho-corrected N and Ciy = uptake * TA
# 2022_sol.27.4_forecast/model.R

# Copyright Iago MOSQUEIRA (WMR), 2022
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(FLasher)

load('data/aap.rda')

# --- SETUP

# TAC advice current year

tac <- FLQuant(9152, dimnames=list(age='all', year=2023), units="tonnes")

# data, intermediate, advice and final YEARS

dy <- dims(run)$maxyear
iy <- dy + 1
ay <- iy + 1
fy <- ay + 1

# GEOMEAN but last 5 years
rec1gm <- exp(mean(log(window(stock.n(run)["1",], end=-5))))

# GENERATE targets from refpts
targets <- expand(as(refpts, 'FLQuant'), year=ay)


# --- SETUP future

# 3 years, 5 years wts/selex, 3 years discards
fut <- stf(run, nyears=3, wts.nyears=3, fbar.nyears=5, disc.nyears=3)

# SSB Mohn's rho correction of stock.n
stock.n(fut)[, '2022'] <- stock.n(run)[, '2022'] * (1 / (1 + mrho$SSB))

# RECALCULATE F
harvest(fut)[,'2022'] <- harvest(stock.n(fut)[,'2022'],
  catch.n(fut)[,'2022'], m(fut)[,'2022']) 

# SET geomean SRR
gmsrr <- predictModel(model=rec~a, params=FLPar(c(rec1gm), units="thousands",
  dimnames=list(params="a", year=seq(iy, length=3), iter=1)))

# --- PROJECT catch options

# Ciy = tac_2023 * mean(official / tac (2020-2022))
# mean(uptake) = 0.469

Ciy <- FLQuants(catch=tac * mean(uptake), fbar=targets["Fmsy",])

fut <- fwd(fut, sr=gmsrr, control=as(Ciy, "fwdControl"))

# DEFINE catch options
catch_options <- list(

  # advice
  advice=FLQuants(fbar=targets["Fmsy",] *
    min((ssb(fut)[, ac(ay)] / refpts$Btrigger), 1)),

  # FMSY
  Fmsy=FLQuants(fbar=targets["Fmsy",]),

  # lowFMSY
  lFmsy=FLQuants(fbar=targets["lFmsy",]),

  # uppFMSY
  uFmsy=FLQuants(fbar=targets["uFmsy",]),

  # lowFMSYadvice
  lFmsyadvice=FLQuants(fbar=targets["lFmsy",] *
    min((ssb(fut)[, ac(ay)] / refpts$Btrigger), 1)),

  # F0
  F0=FLQuants(fbar=FLQuant(0, dimnames=list(age='all', year=ay))),

  # Fpa
  Fpa=FLQuants(fbar=targets["Fpa",]),

  # F05noAR
  F05noAR=FLQuants(fbar=targets["F05noAR",]),

  # Flim
  Flim=FLQuants(fbar=targets["Flim",]),

  # Bpa
  Bpa=FLQuants(ssb_flash=targets["Bpa",]),

  # Blim
  Blim=FLQuants(ssb_flash=targets["Blim",]),

  # MSYBtrigger
  MSYBtrigger=FLQuants(ssb_flash=targets["Btrigger",]),

  # F iy
  lastF=FLQuants(fbar=expand(fbar(fut)[, ac(iy)], year=ay)),
                 
  # Fmp 0.20
  Fmp=FLQuants(fbar=FLQuant(0.20, dimnames=list(age='all', year=ay))),

  # TAC iy
  rotac=FLQuants(catch=expand(tac, year=ay))
)

# C0, fy + 1
F0 <- FLQuants(fbar=FLQuant(0, dimnames=list(age='all', year=fy)))

# CONVERT to fwdControl: intermediate year + catch option + F0

fctls <- lapply(catch_options, function(x)
  as(FLQuants(c(x, F0)), "fwdControl")
)

# RUN!

runs <- FLStocks(lapply(fctls, function(x) fwd(fut, sr=gmsrr, control=x)))

# COMPARE

Map(compare, runs, fctls)

# fit + stock, 3028; stock + fit, 3737

catch(runs$advice)[, '2024']


# SAVE

save(runs, fctls, rec1gm, tac, file="model/runs.rda", compress="xz")
