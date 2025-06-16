import os
import pandas as pd
from core.flowAMOEBA import execFlowAMOEBA
import core.shapefile as shp

# --- User parameters ---
#AREAS1_SHP = "input/Hex37_O.shp"
#AREAS2_SHP = "input/Hex37_D.shp"
#FLOW_FILE   = "input/Flow37xLL.txt"
FLOW_FILE   = "input/Flow10by10.txt"
NEI_LVL    = 12
SIGNIF     = 0.01
MC_REPS    = 500
HEADTAIL   = 0.1

# 1) read spatial polygons
#sf1 = shp.Reader("input/Hex37_O.shp")
sf1 = shp.Reader("input/Ogrid10by10.shp")
shapes1 = sf1.shapes()
AREAS1  = [[shape.points] for shape in shapes1]

#sf2 = shp.Reader("input/Hex37_D.shp")
sf2 = shp.Reader("input/Dgrid10by10.shp")
shapes2 = sf2.shapes()
AREAS2  = [[shape.points] for shape in shapes2]



# 2) read flow file: columns O,D,Value
df = pd.read_csv(FLOW_FILE, sep=r"\s+")
FlowValue = dict(
    zip(
        zip(df["O"].astype(int), df["D"].astype(int)),
        df["Flow"].astype(float),
    )
)


# 3) run flowAMOEBA
result = execFlowAMOEBA(
    AREAS1, AREAS2, FlowValue, NeiLvl=12,
    significance=0.01,
    headtail_pct=0.1,
    mc_reps=500,
)

# 4) write output
out_fn = os.path.join("result","flowAMOEBA_Flow10by10_NeiLvl12_Rook.csv")
os.makedirs(os.path.dirname(out_fn), exist_ok=True)
with open(out_fn, "w") as f:
    f.write(result)
print("Results saved as", out_fn)
