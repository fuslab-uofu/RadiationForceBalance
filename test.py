import RadiationForce as rf
import os.path
import pickle
import comparison_plotter as cp


def goat1(filename):    #harvesting was difficult
    if os.path.exists(filename):
        with open(filename, 'rb') as inF:
            goatSkin = pickle.load(inF)
    else:
        goatSkin = rf.RadiationForce("2023-05-10_SkinRFB/", 3, ["abdomen", "skin_CS", "skin_N"], ["no_sample", "water"], [0.6, 0.4, 0.4])
        goatSkin.change_samplename('skin_CS', 'pre_neck_CS')
        rf.save_object(goatSkin, filename)                         #  CS         alive    perfused                                  
    goatSkin.plot() #skin_N_20pct_T3_39p1W
    goatSkin.comparePlot()
    
    # goatSkin.plotRaw("skin_N_30pct_T2_89p5W")
    # goatSkin.plotRaw("skin_N_30pct_T1_89p6W")
    goatSkin.plotRaw()
    goatSkin.write()

def goat2(filename):
    if os.path.exists(filename):
        with open(filename, 'rb') as inF:
            goatSkin = pickle.load(inF)
    else:                                                #naired    perfused       alive
        goatSkin = rf.RadiationForce("RFB 5_24/", 3, ["abdomen", "preskin_CS", "postskin_CS"], ["water"], [.4, .4, .3])
        goatSkin.change_samplename('preskin_CS', 'pre_neck_CS')
        goatSkin.change_samplename('postskin_CS', 'post_neck_CS') 
        rf.save_object(goatSkin, filename)
                 
    goatSkin.plot()     
    # goatSkin.plotRawSingle("preskin_CS_30pct_T3_90p3W")
    goatSkin.comparePlot()
    goatSkin.plotRaw()
    goatSkin.write()



def goat3(filename):    #harvesting was difficult
    if os.path.exists(filename):
        with open(filename, 'rb') as inF:
            goatSkin = pickle.load(inF)
    else:
        goatSkin = rf.RadiationForce("RFB_6_7/", 3, ["abdomen", "pre_neck_CS", "post_neck_CS"], ["water"], [0.3, 0.4, 0.4])
        rf.save_object(goatSkin, filename)
    s = goatSkin.stats(["pre_neck_CS", "post_neck_CS"])   
    print(s)                                                
    goatSkin.plot() 
    goatSkin.comparePlot()
    # goatSkin.plotRaw("skin_N_30pct_T2_89p5W")
    # goatSkin.plotRaw("skin_N_30pct_T1_89p6W")
    # goatSkin.plotRaw()
    goatSkin.write()

def goat4(filename):    #harvesting was difficult
    if os.path.exists(filename):
        with open(filename, 'rb') as inF:
            goatSkin = pickle.load(inF)
    else:
        goatSkin = rf.RadiationForce("RFB_6_14/", 3, ["abdomen", "pre_neck_CS", "post_neck_CS"], ["water"], [0.3, 0.35, 0.25])
        rf.save_object(goatSkin, filename) 

    s = goatSkin.stats(["pre_neck_CS", "post_neck_CS"], True)
    goatSkin.plot() 
    goatSkin.comparePlot()
    # goatSkin.plotRaw("skin_N_30pct_T2_89p5W")
    # goatSkin.plotRaw("skin_N_30pct_T1_89p6W")
    # goatSkin.plotRaw()
    goatSkin.write()

# goat1("goat1.pkl")
# goat2("goat2.pkl")
# goat3("goat3.pkl")
# goat4("goat4.pkl")

def compare():
    # cp.comparison_plotter(['goat2.pkl', 'goat1.pkl', 'goat3.pkl', 'goat4.pkl'], 'abdomen', color_key=['r', 'g', 'b', 'orange'])
    # cp.comparison_plotter(['goat2.pkl', 'goat1.pkl', 'goat3.pkl', 'goat4.pkl'], 'pre_neck_CS', color_key=['r', 'g', 'b', 'orange'])
    cp.comparison_plotter(['goat2.pkl', 'goat3.pkl', 'goat4.pkl'], 'post_neck_CS', color_key=['r', 'b', 'orange'])

# compare()
    
def agar(filename):
    if os.path.exists(filename):
        with open(filename, 'rb') as inF:
            agar = pickle.load(inF)
    else:
        agar = rf.RadiationForce("RFB 8-7-23/", 3, ["5g", "8g", "10g"], ["Water"], [3, 3, 3])
        rf.save_object(agar, filename)
    # agar.plotRawSingle("5g_20pct_T2_39p5W")
    # agar.plotRawSingle("Water_20pct_T1_39p8W")
    agar.plot(0, 0.5)
    agar.write()

agar("agar.pkl")












#####################################################################################
#cheatsheet


"""
def pig1(filename):
    if os.path.exists(filename):
        with open(filename, 'rb') as inF:
            pigSkin = pickle.load(inF)
    else:
        pigSkin = rf.RadiationForce("RFB_5_26pig/", 3, ["blood"], ["water"], [.3])
        pigSkin.change_samplename('blood', 'neck_hairy')
        rf.save_object(pigSkin, filename)
    pigSkin.plot(0.0, 0.3)
    pigSkin.comparePlot(0.0, 0.3)
    pigSkin.plotRaw(0.0, 0.3)
    pigSkin.write()

"""