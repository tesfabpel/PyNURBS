# @Date:   2016-12-08T20:59:59+00:00
# @Last modified time: 2016-12-08T21:06:13+00:00



dependencies = \
'''This module requires:
     - Numeric Python
'''

import math, time
from Nurbs import Crv as NURBSCrv
try:
    import numpy as np
except ImportError, value:
    print dependencies
    raise

#------------------------------------------------------------------------------#
#              MAIN ENTRY POINT FOR TESTING                                    #
#------------------------------------------------------------------------------#
if __name__ == '__main__':
    cntrl = [[-50., -75., 25., 0., -25.,  75., 50.],
             [25. ,  50., 50., 0., -50., -50., 25.]]
    knots = [0., 0., 0., .2, .4, .6, .8, 1., 1., 1.]
    c = NURBSCrv.Crv(cntrl, knots)
    #c.plot()

#...Define parameter range
    Num = 1e7
    u = np.linspace(0.0, 1.0, num=Num)
    ellapsedTm_alg01 = []
    ellapsedTm_alg02 = []
    speedup_array = []
    error = []


    print 9*"-" + "+" + 36*"-" + "+" + 36*"-" + "+" + 12*"-" + "+" + 12*"-" + "+" + 13*"-" + "+" + 10*"-"
    print "  ui     | Previous algorithm                 | De Boor algorithm                  | Error      | dt (alg 01)| dt (alg 02) | Speed-up  "
    print 9*"-" + "+" + 36*"-" + "+" + 36*"-" + "+" + 12*"-" + "+" + 12*"-" + "+" + 13*"-" + "+" + 10*"-"

    total_starttmr = time.clock()

    for ui in u:
        start_time = time.clock()
        P01 = c.pnt4D([ui])
        dt = (time.clock() - start_time)

        start_time = time.clock()
        P02 = c.pnt4D_deboor([ui])
        dt_deboor = (time.clock() - start_time)
    #...add boath values of ellapsed time to array
        ellapsedTm_alg01.append(dt)
        ellapsedTm_alg02.append(dt_deboor)
        speedup_array.append(dt/dt_deboor)
    #...repack
        P01 = np.array([P01[0][0],P01[1][0],P01[2][0]])
        P02 = np.array([P02[0][0],P02[1][0],P02[2][0]])
        dist = np.linalg.norm(P02-P01)
        error.append(dist)
        print(" %5.5f | %10.5f, %10.5f, %10.5f | %10.5f, %10.5f, %10.5f | %10.5f | %10.5f | %10.5f  | %7.2f " \
              %(ui, P01[0],P01[1],P01[2],\
                    P02[0],P02[1],P02[2], dist, dt*(10**6), dt_deboor*(10**6), dt/dt_deboor ))
    print 9*"-" + "+" + 36*"-" + "+" + 36*"-" + "+" + 12*"-" + "+" + 12*"-" + "+" + 13*"-" + "+" + 10*"-"

    tot_ellapsedtime = (time.clock() - total_starttmr)

#...Get average ellappsed time
    av_elpsdTm_alg01 = np.mean(ellapsedTm_alg01)
    av_elpsdTm_alg02 = np.mean(ellapsedTm_alg02)
    std_elpsdTm_alg01 = 100*np.std(ellapsedTm_alg01)/av_elpsdTm_alg01
    std_elpsdTm_alg02 = 100*np.std(ellapsedTm_alg02)/av_elpsdTm_alg02
    speedup = np.mean(speedup_array)

    print "\n"
    print "-----------------------------------------+------------+-------------------"
    print "                                         | Alg. 01    | Alg. 02 (De Boor) "
    print "-----------------------------------------+------------+-------------------"
    print "        Mean ellapsed time in 10^-6 sec. | %10.3f | %10.3f " %(av_elpsdTm_alg01*(10**6), av_elpsdTm_alg02*(10**6))
    print "                Mean std. deviation in %% | %10.3f | %10.3f " %(std_elpsdTm_alg01, std_elpsdTm_alg02)
    print " Mean speed up factor of deboor agorithm | %10.3f | %10.3f " %(1.0, speedup)
    print "-----------------------------------------+------------+-------------------"
    print "Speed-up factor in %% : %10.3f " %((speedup-1.)*100)
    print "Total ellapsed time: %10.3f seconds " %(tot_ellapsedtime)
    print "Max. error: %10.3f " %(np.amax(error))

    """
    print "\n"
    print "------------------------------------------------------------------------"
    print " Test speed-up of the De Boor aglgorithm if all the u parameter values  "
    print " are given directly as array to method                                  "
    print "------------------------------------------------------------------------"
    startTm = time.clock()
    for ui in u:
        c.pnt4D_deboor([ui])
    dt_single = time.clock() - startTm
    startTm = time.clock()
    c.pnt4D_deboor([ui for ui in u])
    dt_array = time.clock() - startTm

    print " Speed-up factor : %10.5f " %(dt_single/dt_array)
    print "------------------------------------------------------------------------"
    """
