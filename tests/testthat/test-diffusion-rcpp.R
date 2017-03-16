#require(testthat)
context("Diffusion pdiffusion and ddiffusion functions (new rcpp versions).")

is_testing <- TRUE  # Set to TRUE for release package (FALSE only works with extra non-supplied R and C files)

pset1 <- list(a     = 0.8,  
              zr    = 0.3, 
              v     = 3.69, 
              t0    = 0.3, 
              d     = 0, 
              szr   = 0.2, 
              sv    = 0.9, 
              st0   = 0.1,
              bound = "upper", precision = 3)

pset2 <- list(a     = seq(0.8, 0.9, by=0.1), 
              zr    = 0.3, 
              v     = 3.69, 
              t0    = seq(0.3,0.5,length=3), 
              d     = 0, 
              szr   = 0.2, 
              sv    = 0.9, 
              st0   = 0.1,
              bound = c("upper","lower"), precision = 3)

pset3 <- list(a     = 1, 
              zr    = 0.3, 
              v     = 3.69, 
              t0    = seq(0.3,0.5,length=3), 
              d     = 0.1, 
              szr   = 0.2, 
              sv    = 0.9, 
              st0   = 0.1,
              bound = c("upper","lower"), precision = 3)

pset4 <- list(a     = 1, 
              zr    = 0.7, 
              v     = 1.69, 
              t0    = seq(0.3,0.5,length=3), 
              d     = 0, 
              szr   = 0.2, 
              sv    = 0.9, 
              st0   = 0.1,
              bound = c("upper","lower"), precision = 3)

pset5 <- list(a     = 2, 
              zr    = 0.7,
              v     = 1.3,
              t0    = 0.2,
              d     = 0,
              szr   = 0.3,
              sv    = 0.4,
              st0   = 0.5,
              bound = "upper", precision = 3)


# Values to calculate d/p for 
x <- seq (0,4,by=0.1)

if (!is_testing) {
  #### This generates datasets from the old (pre-RCpp) version of rtdists (specifically, from 0.6-6)
  #### That C code is not included in the release package.
  
  # UNLOAD EXISTING RTDISTS, LOAD IN OLD VERSION FOR COMPARISON TESTING
  require (devtools)
  devtools::unload(rtdists)
  #install.packages ("rtdists", lib="tests\\temp_testing\\old_rtdists_0.6-6\\")
  require (rtdists, lib.loc="tests\\temp_testing\\old_rtdists_0.6-6")

  orig_d_set1 <- ddiffusion (x, response=pset1$bound, a=pset1$a, v=pset1$v, t0=pset1$t0, z=pset1$zr, d=pset1$d, sz=pset1$szr, sv=pset1$sv, st0=pset1$st0, precision=pset1$precision)
  orig_d_set2 <- ddiffusion (x, response=pset2$bound, a=pset2$a, v=pset2$v, t0=pset2$t0, z=pset2$zr, d=pset2$d, sz=pset2$szr, sv=pset2$sv, st0=pset2$st0, precision=pset2$precision)
  orig_d_set3 <- ddiffusion (x, response=pset3$bound, a=pset3$a, v=pset3$v, t0=pset3$t0, z=pset3$zr, d=pset3$d, sz=pset3$szr, sv=pset3$sv, st0=pset3$st0, precision=pset3$precision)
  orig_d_set4 <- ddiffusion (x, response=pset4$bound, a=pset4$a, v=pset4$v, t0=pset4$t0, z=pset4$zr, d=pset4$d, sz=pset4$szr, sv=pset4$sv, st0=pset4$st0, precision=pset4$precision)
  orig_d_set5 <- ddiffusion (x, response=pset5$bound, a=pset5$a, v=pset5$v, t0=pset5$t0, z=pset5$zr, d=pset5$d, sz=pset5$szr, sv=pset5$sv, st0=pset5$st0, precision=pset5$precision)
  
  orig_p_set1 <- pdiffusion (x, response=pset1$bound, a=pset1$a, v=pset1$v, t0=pset1$t0, z=pset1$zr, d=pset1$d, sz=pset1$szr, sv=pset1$sv, st0=pset1$st0, precision=pset1$precision)
  orig_p_set2 <- pdiffusion (x, response=pset2$bound, a=pset2$a, v=pset2$v, t0=pset2$t0, z=pset2$zr, d=pset2$d, sz=pset2$szr, sv=pset2$sv, st0=pset2$st0, precision=pset2$precision)
  orig_p_set3 <- pdiffusion (x, response=pset3$bound, a=pset3$a, v=pset3$v, t0=pset3$t0, z=pset3$zr, d=pset3$d, sz=pset3$szr, sv=pset3$sv, st0=pset3$st0, precision=pset3$precision)
  orig_p_set4 <- pdiffusion (x, response=pset4$bound, a=pset4$a, v=pset4$v, t0=pset4$t0, z=pset4$zr, d=pset4$d, sz=pset4$szr, sv=pset4$sv, st0=pset4$st0, precision=pset4$precision)
  orig_p_set5 <- pdiffusion (x, response=pset5$bound, a=pset5$a, v=pset5$v, t0=pset5$t0, z=pset5$zr, d=pset5$d, sz=pset5$szr, sv=pset5$sv, st0=pset5$st0, precision=pset5$precision)
  
  # Print out in a format that makes it easy to just copy/paste below for when is_testing = TRUE  
  tmp <- paste(orig_d_set1, collapse=","); paste0("orig_d_set1 <- c(", tmp, ")")
  tmp <- paste(orig_d_set2, collapse=","); paste0("orig_d_set2 <- c(", tmp, ")")
  tmp <- paste(orig_d_set3, collapse=","); paste0("orig_d_set3 <- c(", tmp, ")")
  tmp <- paste(orig_d_set4, collapse=","); paste0("orig_d_set4 <- c(", tmp, ")")
  tmp <- paste(orig_d_set5, collapse=","); paste0("orig_d_set5 <- c(", tmp, ")")
  
  tmp <- paste(orig_p_set1, collapse=","); paste0("orig_p_set1 <- c(", tmp, ")")
  tmp <- paste(orig_p_set2, collapse=","); paste0("orig_p_set2 <- c(", tmp, ")")
  tmp <- paste(orig_p_set3, collapse=","); paste0("orig_p_set3 <- c(", tmp, ")")
  tmp <- paste(orig_p_set4, collapse=","); paste0("orig_p_set4 <- c(", tmp, ")")
  tmp <- paste(orig_p_set5, collapse=","); paste0("orig_p_set5 <- c(", tmp, ")")
} else {
  orig_d_set1 <- c(0,0,0,0,4.60161372697272,3.041737237641,0.775314844644052,0.20100806320729,0.0554709338440229,0.0161957486935331,0.00496015213247891,0.0015820477670175,0.000522441147958111,0.000177777939582215,6.20919658078304e-05,2.2186660934609e-05,8.08832554025722e-06,3.00144977750857e-06,1.13150079545548e-06,4.32611779274754e-07,1.67506365119956e-07,6.56006229469518e-08,2.59569175886733e-08,1.03669533415991e-08,4.1757727970268e-09,1.69507173542302e-09,6.92977992501777e-10,2.85152550466124e-10,1.18041641819417e-10,4.91349890568927e-11,2.05571224094785e-11,8.64143629954106e-12,3.64849639698318e-12,1.54672293834599e-12,6.58203772075445e-13,2.81091063487899e-13,1.20440317808058e-13,5.17658055663906e-14,2.23139322405947e-14,9.64481885706687e-15,4.17951521546082e-15)
  orig_d_set2 <- c(0,0,0,0,0,0,0.775314844644052,0.0700116011901688,0.775314844644053,0.00353928546940399,0.0161957486935331,0.00353928546940399,0.000522441147958111,0.000251976277630378,0.000522441147958112,2.18490202747983e-05,2.21866609346089e-05,2.18490202747983e-05,1.13150079545548e-06,2.16099430669139e-06,1.13150079545548e-06,2.34429558069607e-07,6.56006229469517e-08,2.34429558069607e-07,4.1757727970268e-09,2.72074791976855e-08,4.17577279702683e-09,3.32247839607158e-09,2.85152550466124e-10,3.32247839607157e-09,2.05571224094785e-11,4.22014357427282e-10,2.05571224094784e-11,5.52966324877327e-11,1.54672293834599e-12,5.52966324877327e-11,1.20440317808058e-13,7.42916764373375e-12,1.20440317808058e-13,1.01875094681317e-12,9.6448188570669e-15)
  orig_d_set3 <- c(0,0,0,0,0.224316377535072,0,0.937541143121596,0.126753167747569,0.937541143121597,0.00757988933744593,0.0455187591671571,0.00757988933744593,0.00330381878930778,0.000731748626620572,0.00330381878930778,8.74395992954433e-05,0.000315516310623682,8.74395992954431e-05,3.62428978864632e-05,1.20187376012792e-05,3.62428978864632e-05,1.82135427331724e-06,4.74235061754887e-06,1.82135427331724e-06,6.82667279031284e-07,2.96284872535752e-07,6.82667279031286e-07,5.08294823925357e-08,1.05615500286041e-07,5.08294823925355e-08,1.72778333632897e-08,9.08472806108125e-09,1.72778333632897e-08,1.67695307663344e-09,2.95409251671028e-09,1.67695307663344e-09,5.23353474367034e-10,3.17669847527765e-10,5.23353474367034e-10,6.14614534542645e-11,9.5452153579115e-11)
  orig_d_set4 <- c(0,0,0,0,0,0,0.855232589686379,0.208607865977777,0.85523258968638,0.0413181064028648,0.120197599664321,0.0413181064028648,0.0200197053245129,0.00804095783640437,0.0200197053245129,0.00161216584359968,0.00357699899429369,0.00161216584359968,0.000669964362901077,0.000330179193913326,0.000669964362901078,6.86827544594374e-05,0.000129768799535354,6.86827544594373e-05,2.57686397845211e-05,1.44568191389061e-05,2.57686397845212e-05,3.07116552815478e-06,5.2151229197864e-06,3.07116552815477e-06,1.07127886505753e-06,6.57262114462798e-07,1.07127886505753e-06,1.41510848321378e-07,2.22696573793223e-07,1.41510848321378e-07,4.67453504638768e-08,3.06204660554486e-08,4.67453504638768e-08,6.65367389909361e-09,9.89119776040409e-09)
  orig_d_set5 <- c(0,0,0,0.000684863555753859,0.0409663461601378,0.171559555851273,0.354416541503711,0.546415314740592,0.724462358219232,0.84105593055722,0.843845451375756,0.772355383476894,0.672361523379911,0.569176374664438,0.474128641968136,0.391267238673863,0.32116664067043,0.262880021576662,0.214906784039583,0.175653423294687,0.143637296268304,0.117562376720732,0.0963335696115949,0.079043613603511,0.064949625185905,0.0534475151592299,0.0440479731300484,0.0363554332579497,0.0300503288312008,0.0248744446422993,0.020618983886394,0.017114922807483,0.0142252500919618,0.0118387378519219,0.00986494693318661,0.0082302221334705,0.00687447921678601,0.00574862452487928,0.0048124799156655,0.00403311157398015,0.00338348192505018)

  orig_p_set1 <- c(0,0,0,0,0.169654748430154,0.644852220438606,0.811721201050489,0.854057849069507,0.865306050934252,0.868482806734998,0.869428910313148,0.869723626727494,0.869819010517484,0.86985091068276,0.869861886388069,0.869865757198307,0.86986715220185,0.869867664636671,0.869867856083721,0.869867928694559,0.869867956607604,0.869867967468629,0.869867971741068,0.869867973438421,0.869867974118822,0.869867974393807,0.869867974505775,0.869867974551679,0.869867974570617,0.869867974578475,0.869867974581754,0.869867974583128,0.869867974583707,0.869867974583952,0.869867974584056,0.8698679745841,0.869867974584119,0.869867974584127,0.869867974584131,0.869867974584132,0.869867974584133)
  orig_p_set2 <- c(0,0,0,0,0,0,0.811721196913999,0.125574562330438,0.811721196913999,0.131997518069757,0.868482802598511,0.131997518069757,0.869819006380996,0.132363879963647,0.869819006380996,0.132391772625441,0.86986575306182,0.132391772625441,0.869867851947233,0.132394304104058,0.869867851947233,0.132394562734449,0.869867963332141,0.132394562734449,0.869867969982335,0.132394591469189,0.869867969982335,0.132394594864978,0.869867970415191,0.132394594864978,0.869867970445266,0.132394595285525,0.869867970445266,0.13239459533954,0.869867970447464,0.13239459533954,0.869867970447631,0.132394595346681,0.869867970447631,0.132394595347648,0.869867970447645)
  orig_p_set3 <- c(0,0,0,0,0.00178292322825537,0,0.775411087834309,0.120922452250843,0.775411087834309,0.132716135127724,0.861195435552079,0.132716135127724,0.865901441795384,0.133576741913059,0.865901441795384,0.133666526168523,0.866276666320826,0.133666526168523,0.866314922623755,0.133677824411988,0.866314922623755,0.133679435165113,0.866319529210655,0.133679435165113,0.866320153470468,0.133679685884491,0.866320153470468,0.13367972750184,0.866320245781243,0.13367972750184,0.866320260352645,0.133679734753962,0.866320260352645,0.133679736066237,0.866320262787705,0.133679736066237,0.866320263198217,0.133679736310875,0.866320263198217,0.133679736357606,0.866320263283258)
  orig_p_set4 <- c(0,0,0,0,0,0,0.781191530696519,0.0504836200207758,0.781191530696518,0.0819600110121512,0.890312005832964,0.0819600110121512,0.906961293667404,0.0880426402075123,0.906961293667404,0.0892403909461952,0.909812660187627,0.0892403909461952,0.910331667455244,0.0894825802843886,0.910331667455244,0.0895324893578018,0.910430168562784,0.0895324893578018,0.910449434981355,0.0895429198770103,0.910449434981355,0.0895451233632549,0.910453289384567,0.0895451233632549,0.910454074007284,0.0895455928255754,0.910454074007284,0.0895456935333022,0.910454235934917,0.0895456935333022,0.910454269724044,0.0895457152583196,0.910454269724044,0.0895457199671402,0.910454276838814)
  orig_p_set5 <- c(0,0,0,7.08446361178149e-06,0.00138436468812972,0.0113507675678673,0.0374267297023063,0.0825133638296158,0.146245531551868,0.225422581216873,0.310520009922219,0.391726173987651,0.464066183870186,0.526107264997353,0.57818038301509,0.62134283021036,0.656860758949882,0.685970510299194,0.70978061414215,0.729242390730573,0.745152325059471,0.758167663671159,0.768826127576991,0.777565489792073,0.784741228585,0.790641683658525,0.795500713500313,0.799508094326108,0.802817980191424,0.805555746159603,0.807823505956937,0.809704554224829,0.811266941644485,0.812566353290477,0.813648428047359,0.814550629885263,0.815303759724723,0.815933178816323,0.816459800292277,0.816900894157634,0.81727074191634)

  test_that("ensure new RCpp ddiffusion function produces the same result as previous C-only versions:", {
  
    tolerance <- 1e-5
  
    expect_equal(orig_d_set1, ddiffusion (x, response=pset1$bound, a=pset1$a, v=pset1$v, t0=pset1$t0, z=pset1$zr, d=pset1$d, sz=pset1$szr, sv=pset1$sv, st0=pset1$st0, precision=pset1$precision), tolerance=tolerance)
    expect_equal(orig_d_set2, ddiffusion (x, response=pset2$bound, a=pset2$a, v=pset2$v, t0=pset2$t0, z=pset2$zr, d=pset2$d, sz=pset2$szr, sv=pset2$sv, st0=pset2$st0, precision=pset2$precision), tolerance=tolerance)
    expect_equal(orig_d_set3, ddiffusion (x, response=pset3$bound, a=pset3$a, v=pset3$v, t0=pset3$t0, z=pset3$zr, d=pset3$d, sz=pset3$szr, sv=pset3$sv, st0=pset3$st0, precision=pset3$precision), tolerance=tolerance)
    expect_equal(orig_d_set4, ddiffusion (x, response=pset4$bound, a=pset4$a, v=pset4$v, t0=pset4$t0, z=pset4$zr, d=pset4$d, sz=pset4$szr, sv=pset4$sv, st0=pset4$st0, precision=pset4$precision), tolerance=tolerance)
    expect_equal(orig_d_set5, ddiffusion (x, response=pset5$bound, a=pset5$a, v=pset5$v, t0=pset5$t0, z=pset5$zr, d=pset5$d, sz=pset5$szr, sv=pset5$sv, st0=pset5$st0, precision=pset5$precision), tolerance=tolerance)
  })
  
  test_that("ensure new RCpp pdiffusion function produces the same result as previous C-only versions:", {
  
    tolerance <- 1e-3   # !! This still produces 1e-2 errors with the 5th parameter set!
    
    expect_equal(orig_p_set1, pdiffusion (x, response=pset1$bound, a=pset1$a, v=pset1$v, t0=pset1$t0, z=pset1$zr, d=pset1$d, sz=pset1$szr, sv=pset1$sv, st0=pset1$st0, precision=pset1$precision), tolerance=tolerance)
    expect_equal(orig_p_set2, pdiffusion (x, response=pset2$bound, a=pset2$a, v=pset2$v, t0=pset2$t0, z=pset2$zr, d=pset2$d, sz=pset2$szr, sv=pset2$sv, st0=pset2$st0, precision=pset2$precision), tolerance=tolerance)
    expect_equal(orig_p_set3, pdiffusion (x, response=pset3$bound, a=pset3$a, v=pset3$v, t0=pset3$t0, z=pset3$zr, d=pset3$d, sz=pset3$szr, sv=pset3$sv, st0=pset3$st0, precision=pset3$precision), tolerance=tolerance)
    expect_equal(orig_p_set4, pdiffusion (x, response=pset4$bound, a=pset4$a, v=pset4$v, t0=pset4$t0, z=pset4$zr, d=pset4$d, sz=pset4$szr, sv=pset4$sv, st0=pset4$st0, precision=pset4$precision), tolerance=tolerance)
    expect_equal(orig_p_set5, pdiffusion (x, response=pset5$bound, a=pset5$a, v=pset5$v, t0=pset5$t0, z=pset5$zr, d=pset5$d, sz=pset5$szr, sv=pset5$sv, st0=pset5$st0, precision=pset5$precision), tolerance=tolerance)
  })

}
