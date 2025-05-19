import sys
import os
from pathlib import Path

import pyubx2.ubxmessage

path = Path(os.path.dirname(os.path.abspath(__file__)))
home = str(path.parent)

import pyubx2
import pynmeagps
import pandas as pd
import numpy as np

import datetime

"""EDIT THIS IF YOU WANT TO CHANGE THE INPUT FILE"""
UBXFILE = "C:\\Users\\jakub\\gnss\\urban_canyon_gnss\\data\\2024_10_7\\raw\\line22_fromhostivartopohorelec.ubx"


ubxfile = UBXFILE


print("Processing file: " + ubxfile)
stream = open(ubxfile, "rb")
ubx = pyubx2.UBXReader(stream, protfilter= pyubx2.UBX_PROTOCOL)



def time_from_rawx(msg):
    tow = msg.rcvTow
    week = msg.week
    leap = msg.leapS
    return datetime.datetime(1980, 1, 6) + datetime.timedelta(weeks=week, seconds=tow)



msgs = []

# for outstream in [out_gps, out_gal, out_bds]:
#     outstream.write("time,leaps \n")

df = pd.DataFrame(columns=["time", "leaps", "gnssId", "svId", "freqId", "code", "phase", "doppler", "locktime", "snr", "codeStd", "phaseStd", "dopplerStd", "codeValid", "phaseValid"])

rows = []
while True:
    msg_raw, msg_parsed = ubx.read()
    if msg_parsed is None:
        break
    elif msg_parsed.msg_id == b'\x15':
        time = time_from_rawx(msg_parsed)
        leaps = msg_parsed.leapS

        for i in range(1,msg_parsed.numMeas+1):
            if i < 10:
                zero = "0"
            else:
                zero = ""
            code = getattr(msg_parsed,'prMes_' + zero + str(i))
            phase = getattr(msg_parsed, 'cpMes_' + zero + str(i))
            doppler = getattr(msg_parsed, 'doMes_' + zero + str(i))
            locktime = getattr(msg_parsed, 'locktime_' + zero + str(i))
            snr = getattr(msg_parsed, 'cno_' + zero + str(i))
            gnssId = getattr(msg_parsed, 'gnssId_' + zero + str(i))
            svId = getattr(msg_parsed, 'svId_' + zero + str(i))
            freqId = getattr(msg_parsed, 'sigId_' + zero + str(i))

            codeStd = getattr(msg_parsed, 'prStd_' + zero + str(i))
            phaseStd = getattr(msg_parsed, 'cpStd_' + zero + str(i))
            dopplerStd = getattr(msg_parsed, 'doStd_' + zero + str(i))
            codeValid = getattr(msg_parsed, 'prValid_' + zero + str(i))
            phaseValid = getattr(msg_parsed, 'cpValid_' + zero + str(i))

            row = {"time": time, "leaps": leaps, "gnssId": gnssId, "svId": svId, "freqId": freqId, "code": code, "phase": phase, "doppler": doppler, "locktime": locktime, "snr": snr, "codeStd": codeStd, "phaseStd": phaseStd, "dopplerStd": dopplerStd, "codeValid": codeValid, "phaseValid": phaseValid}
            rows.append(row)
            # df = pd.concat([df, pd.DataFrame(row, index=[0])], ignore_index=True)

df = pd.DataFrame(rows)
df.to_csv(ubxfile.replace(".ubx", ".csv"), index=False)
          
#         print(dir(msg_parsed))
#         print(msg_parsed._ubxID)
# print(msgs)       

# output.close()
stream.close()




times = df["time"].unique()


gps = []
gal = []
bds = []


tmp = df
gpss = tmp[tmp["gnssId"] == 0]
gals = tmp[tmp["gnssId"] == 2]
bdss = tmp[tmp["gnssId"] == 3]




for time in times:
    print(time)
    tmp = df[df["time"] == time]
    gpss = tmp[tmp["gnssId"] == 0]
    gals = tmp[tmp["gnssId"] == 2]
    bdss = tmp[tmp["gnssId"] == 3]

    for id in bdss.svId.unique():
        leaps = bdss[bdss["svId"] == id].leaps.values[0]
        if 0 in bdss[bdss["svId"] == id].freqId.values:
            C1I = bdss[(bdss["svId"] == id) & (bdss["freqId"] == 0)].code.values[0]
            L1I = bdss[(bdss["svId"] == id) & (bdss["freqId"] == 0)].phase.values[0]
            D1I = bdss[(bdss["svId"] == id) & (bdss["freqId"] == 0)].doppler.values[0]
            C1Istd = bdss[(bdss["svId"] == id) & (bdss["freqId"] == 0)].codeStd.values[0]
            L1Istd = bdss[(bdss["svId"] == id) & (bdss["freqId"] == 0)].phaseStd.values[0]
            D1Istd = bdss[(bdss["svId"] == id) & (bdss["freqId"] == 0)].dopplerStd.values[0]
            C1Ivalid = bdss[(bdss["svId"] == id) & (bdss["freqId"] == 0)].codeValid.values[0]
            L1Ivalid = bdss[(bdss["svId"] == id) & (bdss["freqId"] == 0)].phaseValid.values[0]
            C1Isnr = bdss[(bdss["svId"] == id) & (bdss["freqId"] == 0)].snr.values[0]
        else:
            C1I = np.nan
            L1I = np.nan
            D1I = np.nan
            C1Istd = np.nan
            L1Istd = np.nan
            D1Istd = np.nan
            C1Ivalid = np.nan
            L1Ivalid = np.nan
            C1Isnr = np.nan

        if 2 in bdss[(bdss["svId"] == id)].freqId.values:
            C7I = bdss[(bdss["svId"] == id) & (bdss["freqId"] == 2)].code.values[0]
            L7I = bdss[(bdss["svId"] == id) & (bdss["freqId"] == 2)].phase.values[0]
            D7I = bdss[(bdss["svId"] == id) & (bdss["freqId"] == 2)].doppler.values[0]
            C7Istd = bdss[(bdss["svId"] == id) & (bdss["freqId"] == 2)].codeStd.values[0]
            L7Istd = bdss[(bdss["svId"] == id) & (bdss["freqId"] == 2)].phaseStd.values[0]
            D7Istd = bdss[(bdss["svId"] == id) & (bdss["freqId"] == 2)].dopplerStd.values[0]
            C7Ivalid = bdss[(bdss["svId"] == id) & (bdss["freqId"] == 2)].codeValid.values[0]
            L7Ivalid = bdss[(bdss["svId"] == id) & (bdss["freqId"] == 2)].phaseValid.values[0]
            C7Isnr = bdss[(bdss["svId"] == id) & (bdss["freqId"] == 2)].snr.values[0]
        else:
            C7I = np.nan
            L7I = np.nan
            D7I = np.nan
            C7Istd = np.nan
            L7Istd = np.nan
            D7Istd = np.nan
            C7Ivalid = np.nan
            L7Ivalid = np.nan
            C7Isnr = np.nan

        row = {"time": time,
                "leaps": leaps,
                "gnssId": 3,
                "svId": id,
                "C1I": C1I,
                "L1I": L1I,
                "D1I": D1I,
                "C1Istd": C1Istd,
                "L1Istd": L1Istd,
                "D1Istd": D1Istd,
                "C1Ivalid": C1Ivalid,
                "L1Ivalid": L1Ivalid,
                "C1Isnr": C1Isnr,
                "C7I": C7I,
                "L7I": L7I,
                "D7I": D7I,
                "C7Istd": C7Istd,
                "L7Istd": L7Istd,
                "D7Istd": D7Istd,
                "C7Ivalid": C7Ivalid,
                "L7Ivalid": L7Ivalid,
                "C7Isnr": C7Isnr
                }
        bds.append(row)

    for id in gals.svId.unique():
        leaps = gals[gals["svId"] == id].leaps.values[0]
        if 0 in gals[gals["svId"] == id].freqId.values:
            C1X = gals[(gals["svId"] == id) & (gals["freqId"] == 0)].code.values[0]
            L1X = gals[(gals["svId"] == id) & (gals["freqId"] == 0)].phase.values[0]
            D1X = gals[(gals["svId"] == id) & (gals["freqId"] == 0)].doppler.values[0]
            C1Xstd = gals[(gals["svId"] == id) & (gals["freqId"] == 0)].codeStd.values[0]
            L1Xstd = gals[(gals["svId"] == id) & (gals["freqId"] == 0)].phaseStd.values[0]
            D1Xstd = gals[(gals["svId"] == id) & (gals["freqId"] == 0)].dopplerStd.values[0]
            C1Xvalid = gals[(gals["svId"] == id) & (gals["freqId"] == 0)].codeValid.values[0]
            L1Xvalid = gals[(gals["svId"] == id) & (gals["freqId"] == 0)].phaseValid.values[0]
            C1Xsnr = gals[(gals["svId"] == id) & (gals["freqId"] == 0)].snr.values[0]
        else:
            C1X = np.nan
            L1X = np.nan
            D1X = np.nan
            C1Xstd = np.nan
            L1Xstd = np.nan
            D1Xstd = np.nan
            C1Xvalid = np.nan
            L1Xvalid = np.nan
            C1Xsnr = np.nan

        if 6 in gals[(gals["svId"] == id)].freqId.values:
            C7X = gals[(gals["svId"] == id) & (gals["freqId"] == 6)].code.values[0]
            L7X = gals[(gals["svId"] == id) & (gals["freqId"] == 6)].phase.values[0]
            D7X = gals[(gals["svId"] == id) & (gals["freqId"] == 6)].doppler.values[0]
            C7Xstd = gals[(gals["svId"] == id) & (gals["freqId"] == 6)].codeStd.values[0]
            L7Xstd = gals[(gals["svId"] == id) & (gals["freqId"] == 6)].phaseStd.values[0]
            D7Xstd = gals[(gals["svId"] == id) & (gals["freqId"] == 6)].dopplerStd.values[0]
            C7Xvalid = gals[(gals["svId"] == id) & (gals["freqId"] == 6)].codeValid.values[0]
            L7Xvalid = gals[(gals["svId"] == id) & (gals["freqId"] == 6)].phaseValid.values[0]
            C7Xsnr = gals[(gals["svId"] == id) & (gals["freqId"] == 6)].snr.values[0]
        else:
            C7X = np.nan
            L7X = np.nan
            D7X = np.nan
            C7Xstd = np.nan
            L7Xstd = np.nan
            D7Xstd = np.nan
            C7Xvalid = np.nan
            L7Xvalid = np.nan
            C7Xsnr = np.nan
        row = {"time": time,
                "leaps": leaps,
                "gnssId": 2,
                "svId": id,
                "C1X": C1X,
                "L1X": L1X,
                "D1X": D1X,
                "C1Xstd": C1Xstd,
                "L1Xstd": L1Xstd,
                "D1Xstd": D1Xstd,
                "C1Xvalid": C1Xvalid,
                "L1Xvalid": L1Xvalid,
                "C1Xsnr": C1Xsnr,
                "C7X": C7X,
                "L7X": L7X,
                "D7X": D7X,
                "C7Xstd": C7Xstd,
                "L7Xstd": L7Xstd,
                "D7Xstd": D7Xstd,
                "C7Xvalid": C7Xvalid,
                "L7Xvalid": L7Xvalid,
                "C7Xsnr": C7Xsnr
                }
        gal.append(row)


    for id in gpss.svId.unique():
        leaps = gpss[gpss["svId"] == id].leaps.values[0]
        if 0 in gpss[gpss["svId"] == id].freqId.values:
            C1C = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 0)].code.values[0]
            L1C = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 0)].phase.values[0]
            D1C = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 0)].doppler.values[0]
            C1Cstd = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 0)].codeStd.values[0]
            L1Cstd = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 0)].phaseStd.values[0]
            D1Cstd = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 0)].dopplerStd.values[0]
            C1Cvalid = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 0)].codeValid.values[0]
            L1Cvalid = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 0)].phaseValid.values[0]
            C1Csnr = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 0)].snr.values[0]
        else:
            C1C = np.nan
            L1C = np.nan
            D1C = np.nan
            C1Cstd = np.nan
            L1Cstd = np.nan
            D1Cstd = np.nan
            C1Cvalid = np.nan
            L1Cvalid = np.nan
            C1Csnr = np.nan

        if 3 in gpss[(gpss["svId"] == id)].freqId.values:
            C2X = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 3)].code.values[0]
            L2X = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 3)].phase.values[0]
            D2X = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 3)].doppler.values[0]
            C2Xstd = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 3)].codeStd.values[0]
            L2Xstd = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 3)].phaseStd.values[0]
            D2Xstd = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 3)].dopplerStd.values[0]
            C2Xvalid = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 3)].codeValid.values[0]
            L2Xvalid = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 3)].phaseValid.values[0]
            C2Xsnr = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 3)].snr.values[0]
        elif 4 in gpss[(gpss["svId"] == id)].freqId.values:
            C2X = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 4)].code.values[0]
            L2X = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 4)].phase.values[0]
            D2X = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 4)].doppler.values[0]
            C2Xstd = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 4)].codeStd.values[0]
            L2Xstd = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 4)].phaseStd.values[0]
            D2Xstd = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 4)].dopplerStd.values[0]
            C2Xvalid = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 4)].codeValid.values[0]
            L2Xvalid = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 4)].phaseValid.values[0]
            C2Xsnr = gpss[(gpss["svId"] == id) & (gpss["freqId"] == 4)].snr.values[0]
        else:
            C2X = np.nan
            L2X = np.nan
            D2X = np.nan
            C2Xstd = np.nan
            L2Xstd = np.nan
            D2Xstd = np.nan
            C2Xvalid = np.nan
            L2Xvalid = np.nan
            C2Xsnr = np.nan

        row = {"time": time,
               "leaps": leaps,
                "gnssId": 0,
                "svId": id,
                "C1C": C1C,
                "L1C": L1C,
                "D1C": D1C,
                "C1Cstd": C1Cstd,
                "L1Cstd": L1Cstd,
                "D1Cstd": D1Cstd,
                "C1Cvalid": C1Cvalid,
                "L1Cvalid": L1Cvalid,
                "C1Csnr": C1Csnr,
                "C2X": C2X,
                "L2X": L2X,
                "D2X": D2X,
                "C2Xstd": C2Xstd,
                "L2Xstd": L2Xstd,
                "D2Xstd": D2Xstd,
                "C2Xvalid": C2Xvalid,
                "L2Xvalid": L2Xvalid,
                "C2Xsnr": C2Xsnr
                }
        gps.append(row)

        

gps = pd.DataFrame(gps)
gps.rename(columns = {'time':'Time','svId':'SatelliteID'}, inplace=True)
gps.to_csv(ubxfile.replace("ubx", "gps.csv"), index=False)

gal = pd.DataFrame(gal)
gal.rename(columns = {'time':'Time','svId':'SatelliteID'}, inplace=True)
gal.to_csv(ubxfile.replace("ubx", "gal.csv"), index=False)

bds = pd.DataFrame(bds)
bds.rename(columns = {'time':'Time','svId':'SatelliteID'}, inplace=True)
bds.to_csv(ubxfile.replace("ubx", "bds.csv"), index=False)