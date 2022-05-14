# -*- coding: utf-8 -*-
"""
Created on Fri May 13 16:00:55 2022

@author: jmcti
"""
import matplotlib.pyplot as plt
import numpy as np

tarr = np.arange(0,100)*10000
arr1 = np.asarray([0.7386097999988124, 2.471157399999356, 4.206330400000297, 5.947380400000839, 7.679880200001207, 9.412366000000475, 11.150534100001096, 12.900614700000006, 14.668469000000186, 16.402305100000376, 18.13889560000098, 19.875974900001893, 21.616051799999696, 23.354801400000724, 25.093702300000587, 26.832000099999277, 28.565624800001387, 30.296648399998958, 32.02777509999942, 33.75891870000123, 35.48869270000068, 37.22955590000129, 38.960107900002185, 40.70141760000115, 42.44251540000187, 44.19560420000198, 45.928803300001164, 47.66511869999886, 49.395920599999954, 51.1299528999989, 52.8728873, 54.67015890000039, 56.44195509999918, 58.19630469999902, 59.944298499998695, 61.68273989999943, 63.41554349999933, 65.14683739999964, 66.8807779000017, 68.61531280000054, 70.35144639999999, 72.08594359999915, 73.81491989999995, 75.54615199999898, 77.27842030000102, 79.00948659999995, 80.7423447000001, 82.47543560000122, 84.21528849999959, 85.96005380000133, 87.70310660000177, 89.46134820000225, 91.24965150000207, 93.14206439999907, 95.25993120000203, 97.29688879999958, 99.41052570000102, 101.46490520000225, 103.574268700002, 105.66987169999993, 107.70896279999943, 109.74937440000213, 111.7168672000007, 113.83694200000173, 115.79132419999951, 117.92169479999939, 119.9232047000005, 121.92499750000206, 124.03386049999972, 126.18279640000037, 128.4015528000018, 130.40929770000002, 132.14364660000138, 133.8760035000014, 135.6053510999991, 137.33688619999884, 139.06819149999865, 140.79949580000175, 142.54993359999935, 144.33948880000025, 146.07681660000162, 147.81767730000138, 149.55102650000117, 151.2890993000001, 153.03140700000222, 154.77251220000107, 156.50686290000158, 158.24530420000156, 159.9881501000018, 161.73627530000158, 163.47505490000185, 165.21609890000036, 166.96603070000128, 168.70919680000225, 170.4438783000005, 172.18447799999922, 173.926209199999, 175.6761633999995, 177.4200626999991, 179.17381820000082])
arr2 = np.asarray([0.9212253999976383, 8.076425999999628, 15.236837199998263, 22.436277400000108, 30.380756999998994, 37.872978899999, 45.05431609999869, 52.22959369999808, 59.40515980000055, 66.54596149999998, 73.69899459999942, 80.85141739999744, 88.0303316999998, 95.21814090000044, 102.37046539999938, 109.51153250000061, 116.67702309999731, 123.80935680000039, 130.92605249999906, 138.04155669999818, 145.15887909999947, 152.27691839999898, 159.42429079999783, 166.7062464999981, 175.42850179999732, 185.62840610000057, 195.9117370999993, 205.9762974999976, 215.91566839999723, 226.0789799999984, 236.0063706999972, 244.15828589999728, 251.2801084999992, 258.4024711999991, 265.52576339999723, 272.65726150000046, 279.7748458999995, 286.8900920999986, 294.03346260000035, 301.15905439999915, 308.2791518999984, 315.39651439999943, 322.5196765999972, 329.64388549999785, 336.7636046999978, 343.8987980999991, 351.03518260000055, 358.17219679999835, 365.30106739999974, 372.4200448999982, 379.5395967999975, 386.6624269999993, 393.78534279999803, 400.9058134000006, 408.0283963999973, 415.15453290000005, 422.56038999999873, 429.7719404999989, 436.9110596999999, 444.3486178999992, 451.49431569999797, 458.6559653000004, 465.78851609999765, 472.9316082999976, 480.3286122999998, 488.0231953000002, 495.51216600000043, 502.96473139999944, 512.1967654, 520.9566491999976, 529.6184474999973, 538.6641731999989, 547.737728199998, 555.5954963999975, 562.7447670000001, 569.9102375000002, 577.0822604999994, 584.2776059999997, 591.4063303000003, 598.5438902999995, 605.6694805999978, 612.8044539000002, 620.4369894999982, 629.2905973999987, 638.6262394999976, 647.4850804999987, 655.3172467999975, 662.4634388999984, 669.6001235000003, 676.7359286999999, 683.8987309999975, 691.0413673000003, 698.1687724999974, 705.3776668999999, 712.5726525999999, 719.8177236999982, 727.0459327999997, 734.2530401000004, 741.4205767999993, 748.6021933999982])
plt.plot(tarr,arr1,label='Forward Euler')
plt.plot(tarr,arr2,label='4th Order Runge Kutta')
plt.title('Forward Euler vs. Runge Kutta Simulation Speed')
plt.ylabel('Simulation time')
plt.xlabel('time step')
plt.legend()
plt.savefig('FEvsRKtime.png', dpi=1200)