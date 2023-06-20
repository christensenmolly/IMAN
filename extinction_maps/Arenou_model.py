#!/usr/bin/python
# Arenou model



def _getarenouparams(ll,bb):
  """
  Input: galactic coordinates
  Output: Arenou 1992 alpha, beta, gamma, rr0, saa.
  NOTE: upper limit is taken as <= instead of < 

  @param ll: Galactic Longitude (in degrees)
  @type ll: float
  @param bb: Galactic Lattitude (in degrees)
  @type bb: float
  """
  if -90. <= bb < -60.:
    gamma = 0
    if   0. <= ll < 29.:
      alpha = 2.22534 ; beta = -6.00212 ; rr0 = 0.052 ; saa = 13
    elif 29.  <= ll < 57.:
      alpha = 3.35436 ; beta = -14.74567 ; rr0 = 0.035 ; saa = 40
    elif 57.  <= ll < 85.:
      alpha = 2.77637 ; beta = -9.62706 ; rr0 = 0.042 ; saa = 15
    elif 85.  <= ll < 110.:
      alpha = 4.44190 ; beta = -19.92097 ; rr0 = 0.025 ; saa = 36
    elif 110. <= ll < 150.:
      alpha = 4.46685 ; beta = -26.07305 ; rr0 = 0.026 ; saa = 28
    elif 150. <= ll < 180.:
      alpha = 7.63699 ; beta = -46.10856  ; rr0 = 0.014 ; saa = 38
    elif 180. <= ll < 210.:
      alpha = 2.43412 ; beta = -8.69913 ; rr0 = 0.050 ; saa = 36
    elif 210. <= ll < 240.:
      alpha = 3.34481 ; beta = -13.93228 ; rr0 = 0.035 ; saa = 38
    elif 240. <= ll < 270.:
      alpha = 1.40733 ; beta = -3.43418  ; rr0 = 0.091 ; saa = 30
    elif 270. <= ll < 300.:
      alpha = 1.64466 ; beta = -3.97380 ; rr0 = 0.074 ; saa = 28
    elif 300. <= ll < 330.:
      alpha = 2.12696 ; beta = -6.05682 ; rr0 = 0.056 ; saa = 14
    elif 330. <= ll <= 360.:
      alpha = 2.34636 ; beta = -8.17784 ; rr0 = 0.052 ; saa = 16

  if -60. <= bb < -45.:
    gamma = 0
    if   0.  <= ll < 30.:
      alpha = 2.77060 ; beta = -9.52310 ; rr0 = 0.145 ; saa = 16
    elif 30. <= ll < 60.:
      alpha = 1.96533 ; beta = -5.63251 ; rr0 = 0.174 ; saa = 6
    elif 60. <= ll < 110.:
      alpha = 1.93622 ; beta = -13.31757 ; rr0 = 0.073 ; saa = 26
    elif 110. <= ll < 180.:
      alpha = 1.05414 ; beta = -2.36540 ; rr0 = 0.223 ; saa = 74
    elif 180. <= ll < 210.:
      alpha = 1.39990 ; beta = -1.35325 ; rr0 = 0.252 ; saa = 10
    elif 210. <= ll < 240.:
      alpha = 2.73481 ; beta = -11.70266 ; rr0 = 0.117 ; saa = 8
    elif 240. <= ll < 270.:
      alpha = 2.99784 ; beta = -11.64272 ; rr0 = 0.129 ; saa = 3
    elif 270. <= ll < 300.:
      alpha = 3.23802 ; beta = -11.63810 ; rr0 = 0.139 ; saa = 7
    elif 300. <= ll < 330.:
      alpha = 1.72679 ; beta = -6.05085 ; rr0 = 0.143 ; saa = 7
    elif 330. <= ll <= 360.:
      alpha = 1.88890 ; beta = -5.51861 ; rr0 = 0.171 ; saa = 14

  if -45. <= bb < -30.:
    gamma = 0
    if 0.  <= ll < 30.:
      alpha = 1.98973 ; beta = -4.86206  ; rr0 = 0.205 ; saa = 6
    elif 30.  <= ll < 60.:
      alpha = 1.49901 ; beta = -3.75837  ; rr0 = 0.199 ; saa = 16
    elif 60.  <= ll < 90.:
      alpha = 0.90091 ; beta = -1.30459  ; rr0 = 0.329 ; saa = 73
    elif 90.  <= ll < 120.:
      alpha = 1.94200 ; beta = -6.26833  ; rr0 = 0.155 ; saa = 18
    elif 120. <= ll < 160.:
      alpha = -0.37804 ; beta = 10.77372  ; rr0 = 0.210 ; saa = 100
    elif 160. <= ll < 200.:
      alpha = -0.15710 ; beta = 5.03190   ; rr0 = 0.294 ; saa = 24
    elif 200. <= ll < 235.:
      alpha = 3.20162 ; beta = -10.59297 ; rr0 = 0.151 ; saa = 9
    elif 235. <= ll < 265.:
      alpha = 1.95079 ; beta = -4.73280   ; rr0 = 0.206 ; saa = 21
    elif 265. <= ll < 300.:
      alpha = 1.91200 ; beta = -4.97640   ; rr0 = 0.192 ; saa = 17
    elif 300. <= ll < 330.:
      alpha = 2.50487 ; beta = -8.63106  ; rr0 = 0.145 ; saa = 28
    elif 330. <= ll <= 360.:
      alpha = 2.44394 ; beta = -9.17612  ; rr0 = 0.133 ; saa = 7

  if -30. <= bb < -15.:
    gamma = 0
    if 0. <= ll < 20.:
      alpha = 2.82440 ; beta = -4.78217 ; rr0 = 0.295 ; saa = 32
    elif 20.  <= ll < 40.:
      alpha = 3.84362 ; beta = -8.04690 ; rr0 = 0.239 ; saa = 46
    elif 40.  <= ll < 80.:
      alpha = 0.60365 ; beta = 0.07893  ; rr0 = 0.523 ; saa = 22
    elif 80.  <= ll < 100.:
      alpha = 0.58307 ; beta = -0.21053 ; rr0 = 0.523 ; saa = 53
    elif 100. <= ll < 120.:
      alpha = 2.03861 ; beta = -4.40843 ; rr0 = 0.231 ; saa = 60
    elif 120. <= ll < 140.:
      alpha = 1.14271 ; beta = -1.35635 ; rr0 = 0.421 ; saa = 34
    elif 140. <= ll < 160.:
      alpha = 0.79908 ; beta = 1.48074  ; rr0 = 0.523 ; saa = 61
    elif 160. <= ll < 180.:
      alpha = 0.94260 ; beta = 8.16346  ; rr0 = 0.441 ; saa = 42
    elif 180. <= ll < 200.:
      alpha = 1.66398 ; beta = 0.26775  ; rr0 = 0.523 ; saa = 42
    elif 200. <= ll < 220.:
      alpha = 1.08760 ; beta = -1.02443 ; rr0 = 0.523 ; saa = 45
    elif 220. <= ll < 240.:
      alpha = 1.20087 ; beta = -2.45407 ; rr0 = 0.245 ; saa = 6
    elif 240. <= ll < 260.:
      alpha = 1.13147 ; beta = -1.87916 ; rr0 = 0.301 ; saa = 16
    elif 260. <= ll < 280.:
      alpha = 1.97804 ; beta = -2.92838 ; rr0 = 0.338 ; saa = 21
    elif 280. <= ll < 300.:
      alpha = 1.40086 ; beta = -1.12403 ; rr0 = 0.523 ; saa = 19
    elif 300. <= ll < 320.:
      alpha = 2.06355 ; beta = -3.68278 ; rr0 = 0.280 ; saa = 42
    elif 320. <= ll < 340.:
      alpha = 1.59260 ; beta = -2.18754 ; rr0 = 0.364 ; saa = 15
    elif 340. <= ll <= 360.:
      alpha = 1.45589 ; beta = -1.90598 ; rr0 = 0.382 ; saa = 21

  if -15. <= bb < -5.:
    gamma = 0
    if 0. <= ll < 10.:
      alpha = 2.56438 ; beta = -2.31586 ; rr0 = 0.554 ; saa = 37
    elif 10.  <= ll < 20.:
      alpha = 3.24095 ; beta = -2.78217 ; rr0 = 0.582 ; saa = 38
    elif 20.  <= ll < 30.:
      alpha = 2.95627 ; beta = -2.57422 ; rr0 = 0.574 ; saa = 41 ; gamma = 0.08336
    elif 30.  <= ll < 40.:
      alpha = 1.85158 ; beta = -0.67716 ; rr0 = 1.152 ; saa = 4
    elif 40.  <= ll < 50.:
      alpha = 1.60770 ; beta = 0.35279  ; rr0 = 0.661 ; saa = 24
    elif 50.  <= ll < 60.:
      alpha = 0.69920 ; beta = -0.09146 ; rr0 = 0.952 ; saa = 2  ; gamma = 0.12839
    elif 60.  <= ll < 80.:
      alpha = 1.36189 ; beta = -1.05290 ; rr0 = 0.647 ; saa = 45 ; gamma = 0.16258
    elif 80.  <= ll < 90.:
      alpha = 0.33179 ; beta = 0.37629  ; rr0 = 1.152 ; saa = 62
    elif 90.  <= ll < 100.:
      alpha = 1.70303 ; beta = -0.75246 ; rr0 = 1.132 ; saa = 31
    elif 100. <= ll < 110.:
      alpha = 1.97414 ; beta = -1.59784 ; rr0 = 0.618 ; saa = 35 ; gamma = 0.12847
    elif 110. <= ll < 120.:
      alpha = 1.07407 ; beta = -0.40066 ; rr0 = 1.152 ; saa = 14 ; gamma = 0.17698
    elif 120. <= ll < 130.:
      alpha = 1.69495 ; beta = -1.00071 ; rr0 = 0.847 ; saa = 28 ; gamma = 0.08567
    elif 130. <= ll < 140.:
      alpha = 1.51449 ; beta = -0.08441 ; rr0 = 0.897 ; saa = 12
    elif 140. <= ll < 150.:
      alpha = 1.87894 ; beta = -0.73314 ; rr0 = 1.152 ; saa = 23
    elif 150. <= ll < 160.:
      alpha = 1.43670 ; beta = 0.67706  ; rr0 = 0.778 ; saa = 3  ; gamma = 0.42624
    elif 160. <= ll < 180.:
      alpha = 6.84802 ; beta = -5.06864 ; rr0 = 0.676 ; saa = 50
    elif 180. <= ll < 190.:
      alpha = 4.16321 ; beta = -5.80016 ; rr0 = 0.359 ; saa = 51 ; gamma = 0.60387
    elif 190. <= ll < 200.:
      alpha = 0.78135 ; beta = -0.27826 ; rr0 = 1.152 ; saa = 4
    elif 200. <= ll < 210.:
      alpha = 0.85535 ; beta = 0.20848  ; rr0 = 1.152 ; saa = 17
    elif 210. <= ll < 220.:
      alpha = 0.52521 ; beta = 0.65726  ; rr0 = 1.152 ; saa = 7
    elif 220. <= ll < 230.:
      alpha = 0.88376 ; beta = -0.44519 ; rr0 = 0.993 ; saa = 40 ; gamma = 0.06013
    elif 230. <= ll < 240.:
      alpha = 0.42228 ; beta = -0.26304  ; rr0 = 0.803 ; saa = 26
    elif 240. <= ll < 250.:
      alpha = 0.71318 ; beta = -0.67229 ; rr0 = 0.530 ; saa = 55 ; gamma = 0.20994
    elif 250. <= ll < 260.:
      alpha = 0.99606 ; beta = -0.70103 ; rr0 = 0.710 ; saa = 48 ; gamma = 0.01323
    elif 260. <= ll < 270.:
      alpha = 0.91519 ; beta = -0.39690 ; rr0 = 1.152 ; saa = 48 ; gamma = 0.01961
    elif 270. <= ll < 280.:
      alpha = 0.85791 ; beta = -0.29115 ; rr0 = 1.152 ; saa = 19
    elif 280. <= ll < 290.:
      alpha = 1.44226 ; beta = -1.09775 ; rr0 = 0.657 ; saa = 39
    elif 290. <= ll < 300.:
      alpha = 2.55486 ; beta = -1.68293 ; rr0 = 0.759 ; saa = 31
    elif 300. <= ll < 310.:
      alpha = 3.18047 ; beta = -2.69796 ; rr0 = 0.589 ; saa = 40
    elif 310. <= ll < 320.:
      alpha = 2.11235 ; beta = -1.77506 ; rr0 = 0.595 ; saa = 29
    elif 320. <= ll < 330.:
      alpha = 1.75884 ; beta = -1.38574 ; rr0 = 0.635 ; saa = 25
    elif 330. <= ll < 340.:
      alpha = 1.97257 ; beta = -1.55545 ; rr0 = 0.634 ; saa = 34 ; gamma = 0.00043
    elif 340. <= ll < 350.:
      alpha = 1.41497 ; beta = -1.05722 ; rr0 = 0.669 ; saa = 46 ; gamma = 0.03264
    elif 350. <= ll <= 360.:
      alpha = 1.17795 ; beta = -0.95012 ; rr0 = 0.620 ; saa = 46 ; gamma = 0.03339

  if -5. <= bb < 5.:
    gamma = 0
    if 0.   <= ll < 10.:
      alpha = 2.62556 ; beta = -1.11097 ; rr0 = 1.182 ; saa = 57 ; gamma = 0.00340
    elif 10.  <= ll < 20.:
      alpha = 3.14461 ; beta = -1.01140 ; rr0 = 1.555 ; saa = 42
    elif 20.  <= ll < 30.:
      alpha = 4.26624 ; beta = -1.61242 ; rr0 = 1.323 ; saa = 34
    elif 30.  <= ll < 40.:
      alpha = 2.54447 ; beta = -0.12771 ; rr0 = 1.300 ; saa = 30
    elif 40.  <= ll < 50.:
      alpha = 2.27030 ; beta = -0.68720 ; rr0 = 1.652 ; saa = 52 ; gamma = 0.04928
    elif 50.  <= ll < 60.:
      alpha = 1.34359 ; beta = -0.05416 ; rr0 = 2.000 ; saa = 32
    elif 60.  <= ll < 70.:
      alpha = 1.76327 ; beta = -0.26407 ; rr0 = 2.000 ; saa = 37
    elif 70.  <= ll < 80.:
      alpha = 2.20666 ; beta = -0.41651 ; rr0 = 2.000 ; saa = 36
    elif 80.  <= ll < 90.:
      alpha = 1.50130 ; beta = -0.09855 ; rr0 = 1.475 ; saa = 45
    elif 90.  <= ll < 100.:
      alpha = 2.43965 ; beta = -0.82128 ; rr0 = 1.485 ; saa = 36 ; gamma = 0.01959
    elif 100. <= ll < 110.:
      alpha = 3.35775 ; beta = -1.16400 ; rr0 = 0.841 ; saa = 35 ; gamma = 0.00298
    elif 110. <= ll < 120.:
      alpha = 2.60621 ; beta = -0.68687 ; rr0 = 1.897 ; saa = 36
    elif 120. <= ll < 130.:
      alpha = 2.90112 ; beta = -0.97988 ; rr0 = 1.480 ; saa = 32
    elif 130. <= ll < 140.:
      alpha = 2.55377 ; beta = -0.71214 ; rr0 = 1.793 ; saa = 38
    elif 140. <= ll < 150.:
      alpha = 3.12598 ; beta = -1.21437 ; rr0 = 1.287 ; saa = 23 ; gamma = 0.15298
    elif 150. <= ll < 160.:
      alpha = 3.66930 ; beta = -2.29731 ; rr0 = 0.799 ; saa = 40 ; gamma = 0.33473
    elif 160. <= ll < 170.:
      alpha = 2.15465 ; beta = -0.70690 ; rr0 = 1.524 ; saa = 37 ; gamma = 0.14017
    elif 170. <= ll < 180.:
      alpha = 1.82465 ; beta = -0.60223 ; rr0 = 1.515 ; saa = 29 ; gamma = 0.20730
    elif 180. <= ll < 190.:
      alpha = 1.76269 ; beta = -0.35945 ; rr0 = 2.000 ; saa = 28 ; gamma = 0.08052
    elif 190. <= ll < 200.:
      alpha = 1.06085 ; beta = -0.14211 ; rr0 = 2.000 ; saa = 48
    elif 200. <= ll < 210.:
      alpha = 1.21333 ; beta = -0.23225 ; rr0 = 2.000 ; saa = 57
    elif 210. <= ll < 220.:
      alpha = 0.58326 ; beta = -0.06097 ; rr0 = 2.000 ; saa = 30 ; gamma = 0.36962
    elif 220. <= ll < 230.:
      alpha = 0.74200 ; beta = -0.19293 ; rr0 = 1.923 ; saa = 48 ; gamma = 0.07459
    elif 230. <= ll < 240.:
      alpha = 0.67520 ; beta = -0.21041 ; rr0 = 1.604 ; saa = 49 ; gamma = 0.16602
    elif 240. <= ll < 250.:
      alpha = 0.62609 ; beta = -0.25312 ; rr0 = 1.237 ; saa = 73 ; gamma = 0.14437
    elif 250. <= ll < 260.:
      alpha = 0.61415 ; beta = -0.13788 ; rr0 = 2.000 ; saa = 42 ; gamma = 0.26859
    elif 260. <= ll < 270.:
      alpha = 0.58108 ; beta = 0.01195  ; rr0 = 2.000 ; saa = 40 ; gamma = 0.07661
    elif 270. <= ll < 280.:
      alpha = 0.68352 ; beta = -0.10743 ; rr0 = 2.000 ; saa = 50 ; gamma = 0.00849
    elif 280. <= ll < 290.:
      alpha = 0.61747 ; beta = 0.02675  ; rr0 = 2.000 ; saa = 49
    elif 290. <= ll < 300.:
      alpha = 1.06827 ; beta = -0.26290 ; rr0 = 2.000 ; saa = 44
    elif 300. <= ll < 310.:
      alpha = 1.53631 ; beta = -0.36833 ; rr0 = 2.000 ; saa = 37 ; gamma = 0.02960
    elif 310. <= ll < 320.:
      alpha = 1.94300 ; beta = -0.71445 ; rr0 = 1.360 ; saa = 36 ; gamma = 0.15643
    elif 320. <= ll < 330.:
      alpha = 1.22185 ; beta = -0.40185 ; rr0 = 1.520 ; saa = 48 ; gamma = 0.07354
    elif 330. <= ll < 340.:
      alpha = 1.79429 ; beta = -0.48657 ; rr0 = 1.844 ; saa = 43
    elif 340. <= ll < 350.:
      alpha = 2.29545 ; beta = -0.84096 ; rr0 = 1.365 ; saa = 32
    elif 350. <= ll <= 360.:
      alpha = 2.07408 ; beta = -0.64745 ; rr0 = 1.602 ; saa = 36 ; gamma = 0.12750


  if 5. <= bb < 15.:
    gamma = 0
    if 0.   <= ll < 10.:
      alpha = 2.94213 ; beta = -2.09158 ; rr0 = 0.703 ; saa = 41 ; gamma = 0.05490
    elif 10.  <= ll < 30.:
      alpha = 3.04627 ; beta = 7.71159  ; rr0 = 0.355 ; saa = 45
    elif 30.  <= ll < 40.:
      alpha = 3.78033 ; beta = -3.91956 ; rr0 = 0.482 ; saa = 42
    elif 40.  <= ll < 50.:
      alpha = 2.18119 ; beta = -2.4050  ; rr0 = 0.453 ; saa = 27
    elif 50.  <= ll < 60.:
      alpha = 1.45372 ; beta = -0.49522 ; rr0 = 1.152 ; saa = 31
    elif 60.  <= ll < 70.:
      alpha = 1.05051 ; beta = -1.01704 ; rr0 = 0.516 ; saa = 2
    elif 70. <= ll < 80.:
      alpha = 0.48416 ; beta = -0.27182 ; rr0 = 0.891 ; saa = 94 ; gamma = 0.08639
    elif 80.  <= ll < 90.:
      alpha = 0.61963 ; beta = 0.41697  ; rr0 = 1.152 ; saa = 35 ; gamma = 0.47171
    elif 90.  <= ll < 100.:
      alpha = 4.40348 ; beta = -2.95611 ; rr0 = 0.745 ; saa = 52
    elif 100. <= ll < 120.:
      alpha = 2.50938 ; beta = -0.56541 ; rr0 = 1.152 ; saa = 27
    elif 120. <= ll < 130.:
      alpha = 0.44180 ; beta = 1.58923  ; rr0 = 0.949 ; saa = 4
    elif 130. <= ll < 140.:
      alpha = 3.96084 ; beta = -3.37374 ; rr0 = 0.587 ; saa = 40 ; gamma = 0.34109
    elif 140. <= ll < 160.:
      alpha = 2.53335 ; beta = -0.40541 ; rr0 = 1.152 ; saa = 38
    elif 160. <= ll < 170.:
      alpha = 2.03760 ; beta = -0.66317 ; rr0 = 1.152 ; saa = 23
    elif 170. <= ll < 200.:
      alpha = 1.06946 ; beta = -0.87395 ; rr0 = 0.612 ; saa = 29 ; gamma = 0.29230
    elif 200. <= ll < 210.:
      alpha = 0.86348 ; beta = -0.65870 ; rr0 = 0.655 ; saa = 79 ; gamma = 0.09089
    elif 210. <= ll < 230.:
      alpha = 0.30117 ; beta = -0.16136 ; rr0 = 0.933 ; saa = 17 ; gamma = 0.07495
    elif 230. <= ll < 240.:
      alpha = 0.75171 ; beta = -0.57143 ; rr0 = 0.658 ; saa = 12 ; gamma = 0.00534
    elif 240. <= ll < 250.:
      alpha = 1.97427 ; beta = -2.02654 ; rr0 = 0.487 ; saa = 67
    elif 250. <= ll < 260.:
      alpha = 1.25208 ; beta = -1.47763 ; rr0 = 0.424 ; saa = 19 ; gamma = 0.31600
    elif 260. <= ll < 270.:
      alpha = 0.89448 ; beta = -0.43870 ; rr0 = 1.019 ; saa = 5
    elif 270. <= ll < 280.:
      alpha = 0.81141 ; beta = -0.51001 ; rr0 = 0.795 ; saa = 27 ; gamma = 0.03505
    elif 280. <= ll < 290.:
      alpha = 0.83781 ; beta = -0.44138 ; rr0 = 0.949 ; saa = 50 ; gamma = 0.02820
    elif 290. <= ll < 300.:
      alpha = 1.10600 ; beta = -0.86263 ; rr0 = 0.641 ; saa = 28 ; gamma = 0.03402
    elif 300. <= ll < 310.:
      alpha = 1.37040 ; beta = -1.02779 ; rr0 = 0.667 ; saa = 28 ; gamma = 0.05608
    elif 310. <= ll < 320.:
      alpha = 1.77590 ; beta = -1.26951 ; rr0 = 0.699 ; saa = 37 ; gamma = 0.06972
    elif 320. <= ll < 330.:
      alpha = 1.20865 ; beta = -0.70679 ; rr0 = 0.855 ; saa = 35 ; gamma = 0.02902
    elif 330. <= ll < 340.:
      alpha = 2.28830 ; beta = -1.71890 ; rr0 = 0.666 ; saa = 42 ; gamma = 0.22887
    elif 340. <= ll < 350.:
      alpha = 3.26278 ; beta = -0.94181 ; rr0 = 1.152 ; saa = 38
    elif 350. <= ll <= 360.:
      alpha = 2.58100 ; beta = -1.69237 ; rr0 = 0.763 ; saa = 53

  if 15. <= bb < 30.:
    gamma = 0
    if 0.   <= ll < 20.:
      alpha = 6.23279  ; beta = -10.30384 ; rr0 = 0.302 ; saa = 42
    elif 20.  <= ll < 40.:
      alpha = 4.47693 ; beta = -7.28366  ; rr0 = 0.307 ; saa = 29
    elif 40.  <= ll < 60. :
      alpha =  1.22938 ; beta = -1.19030  ; rr0 = 0.516 ; saa = 5
    elif 60.  <= ll < 80. :
      alpha =  0.84291 ; beta = -1.59338  ; rr0 = 0.265 ; saa = 4
    elif 80.  <= ll < 100. :
      alpha =  0.23996 ; beta = 0.06304   ; rr0 = 0.523 ; saa = 32
    elif 100. <= ll < 140. :
      alpha =  0.40062 ; beta = -1.75628  ; rr0 = 0.114 ; saa = 16
    elif 140. <= ll < 180. :
      alpha =  0.56898 ; beta = -0.53331  ; rr0 = 0.523 ; saa = 41
    elif 180. <= ll < 200. :
      alpha = -0.95721 ; beta = 11.69217   ; rr0 = 0.240 ; saa = 2
    elif 200. <= ll < 220. :
      alpha = -0.19051 ; beta = 1.45670   ; rr0 = 0.376 ; saa = 1
    elif 220. <= ll < 240. :
      alpha =  2.31305 ; beta = -7.82531  ; rr0 = 0.148 ; saa = 95
    elif 240. <= ll < 260.:
      alpha =  1.39169 ; beta = -1.72984  ; rr0 = 0.402 ; saa = 6
    elif 260. <= ll < 280.:
      alpha =  1.59418 ; beta = -1.28296  ; rr0 = 0.523 ; saa = 36
    elif 280. <= ll < 300. :
      alpha =  1.57082 ; beta = -1.97295   ; rr0 = 0.398 ; saa = 10
    elif 300. <= ll < 320. :
      alpha =  1.95998 ; beta = -3.26159  ; rr0 = 0.300 ; saa = 11
    elif 320. <= ll < 340.:
      alpha =  2.59567 ; beta = -4.84133  ; rr0 = 0.268 ; saa = 37
    elif 340. <= ll <= 360.:
      alpha =  5.30273 ; beta = -7.43033  ; rr0 = 0.357 ; saa = 37

  if 30. <= bb < 45.:
    gamma = 0
    if 0.   <= ll < 20.:
      alpha =  2.93960 ; beta  = -6.48049  ; rr0 = 0.227 ; saa = 77
    elif 20.  <= ll < 50.:
      alpha =  1.65864 ; beta  = -9.99317  ; rr0 = 0.083 ; saa = 99
    elif 50.  <= ll < 80.:
      alpha =  1.71831 ; beta  = -7.25286  ; rr0 = 0.118 ; saa = 28
    elif 80.  <= ll < 110.:
      alpha =  1.33617 ; beta  = -10.39799 ; rr0 = 0.064 ; saa = 99
    elif 110. <= ll < 160.:
      alpha = -0.31330 ; beta  = 1.35622   ; rr0 = 0.329 ; saa = 24
    elif 160. <= ll < 190.:
      alpha =  1.51984 ; beta  = -8.69502  ; rr0 = 0.087 ; saa = 99
    elif 190. <= ll < 220.:
      alpha = -0.50758 ; beta  = 4.73320   ; rr0 = 0.250 ; saa = 78
    elif 220. <= ll < 250.:
      alpha =  1.25864 ; beta  = -12.59627 ; rr0 = 0.050 ; saa = 70
    elif 250. <= ll < 280.:
      alpha =  1.54243 ; beta  = -3.76065  ; rr0 = 0.205 ; saa = 10
    elif 280. <= ll < 320.:
      alpha =  2.72258 ; beta  = -7.47806  ; rr0 = 0.182 ; saa = 5
    elif 320. <= ll < 340.:
      alpha =  2.81545 ; beta  = -5.52139  ; rr0 = 0.255 ; saa = 10
    elif 340. <= ll <= 360.:
      alpha =  2.23818 ; beta  = 0.81772   ; rr0 = 0.329 ; saa = 19

  if 45. <= bb < 60.:
    gamma = 0
    if 0.   <= ll < 60.:
      alpha = 1.38587 ; beta  = -9.06536  ; rr0 = 0.076 ; saa = 3
    elif 60.  <= ll < 90.:
      alpha = 2.28570 ; beta  = -9.88812  ; rr0 = 0.116 ; saa = 3
    elif 90.  <= ll < 110.:
      alpha = 1.36385 ; beta  = -8.10127  ; rr0 = 0.084 ; saa = 4
    elif 110. <= ll < 170.:
      alpha = 0.05943 ; beta  = -1.08126  ; rr0 = 0.027 ; saa = 50
    elif 170. <= ll < 200.:
      alpha = 1.40171 ; beta  = -3.21783  ; rr0 = 0.218 ; saa = 99
    elif 200. <= ll < 230.:
      alpha = 0.14718 ; beta  = 3.92670   ; rr0 = 0.252 ; saa = 14
    elif 230. <= ll < 290.:
      alpha = 0.57124 ; beta  = -4.30242  ; rr0 = 0.066 ; saa = 10
    elif 290. <= ll < 330.:
      alpha = 3.69891 ; beta  = -19.62204 ; rr0 = 0.094 ; saa = 5
    elif 330. <= ll <= 360.:
      alpha = 1.19568 ; beta  = -0.45043  ; rr0 = 0.252 ; saa = 9

  if 60. <= bb <= 90.:
    gamma = 0
    if 0.   <= ll < 30.:
      alpha = 0.69443 ;  beta = -0.27600  ; rr0 = 0.153 ; saa = 99
    elif 30.  <= ll < 60.:
      alpha = 1.11811 ;  beta = 0.71179   ; rr0 = 0.085 ; saa = 73
    elif 60.  <= ll < 90.:
      alpha = 1.10427 ;  beta = -2.37654  ; rr0 = 0.123 ; saa = 99
    elif 90.  <= ll < 120.:
      alpha = -0.42211 ; beta = 5.24037   ; rr0 = 0.184 ; saa = 12
    elif 120. <= ll < 150.:
      alpha = 0.87576 ;  beta = -4.38033  ; rr0 = 0.100 ; saa = 35
    elif 150. <= ll < 180.:
      alpha = 1.27477 ;  beta = -4.98307  ; rr0 = 0.128 ; saa = 72
    elif 180. <= ll < 210.:
      alpha = 1.19512 ;  beta = -6.58464  ; rr0 = 0.091 ; saa = 49
    elif 210. <= ll < 240.:
      alpha = 0.97581 ;  beta = -4.89869  ; rr0 = 0.100 ; saa = 95
    elif 240. <= ll < 270.:
      alpha = 0.54379 ;  beta = -0.84403  ; rr0 = 0.207 ; saa = 35
    elif 270. <= ll < 300.:
      alpha = -0.85054 ;  beta = 13.01249  ; rr0 = 0.126 ; saa = 39
    elif 300. <= ll < 330.:
      alpha = 0.74347 ;  beta = -1.39825   ; rr0 = 0.207 ; saa = 10
    elif 330. <= ll <= 360.:
      alpha = 0.77310 ;  beta = -4.45005  ; rr0 = 0.087 ; saa = 16
  #print alpha, beta, gamma, rr0, saa
  return alpha, beta, gamma, rr0, saa



def Arenou(pts):
  """
  Find the predicted V-band extinction (Av) according to the
  3D model for galactic extinction of Arenou et al, "Modelling
  the Galactic interstellar extinction distribution in three
  dimensions", Arenou et al, "A tridimensional model of the
  galactic interstellar extinction" published in Astronomy and
  Astrophysics (ISSN 0004-6361), vol. 258, no. 1, p. 104-111, 1992.
  If the distance is not given, we use the maximal distance r0 for
  that line of sight, as given in the Appendix of Arenou et al. (1992).
  Distance should be given in pc!
  Taken from IVS, KuLeuven: 
  ftp://ftp.ster.kuleuven.be/dist/ehsan/ivs-repos-for-outsiders/ivs/doc/html/ivs.sed.extinctionmodels-pysrc.html#findext_drimmel
  """
  result = []
  for k in range(len(pts)):
      #print k
      [ll,bb,distance] = pts[k]

      # find the Arenou paramaters in the Appendix of Arenou et al. (1992)
      alpha, beta, gamma, rr0, saa = _getarenouparams(ll, bb)


      # compute the visual extinction from the Arenou paramaters using Equation 5
      # and 5bis
      if distance is None:
        Av = alpha*rr0 + beta*rr0**2.
      else:
        distance = distance/1e3 # to kpc
        if distance <= rr0:
          Av = alpha*distance + beta*distance**2.
        else:
          #print alpha,ll,bb
          Av = alpha*rr0 + beta*rr0**2. + (distance-rr0)*gamma
      
      if Av<0.:
          Av=0.
      result.append( [round(Av,3), round(Av*saa/100.,3)])
  return result

#Av = Arenou([[192.38,-42.92,1.609]])
#print(Av)
#exit()