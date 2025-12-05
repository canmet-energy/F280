import pandas as pd
import numpy as np
import os.path
import pathlib
import math

def getFoundationWeather():
    sPathDbs = str(pathlib.Path(__file__).resolve().parent).__str__()
    sPathDbs=sPathDbs.replace('\\','/')
    sWthFile = sPathDbs+'/FdnWth.csv'
    df = pd.read_csv(sWthFile, encoding = "ISO-8859-1")
    return df
    
class basesimp:
    def __init__(self,iLoc=187,config='BCIN_4',length=12.4,width=6.4,overlap=0.6,rsi_interior=0.0, rsi_exterior=0.0, rsi_slab=0.0, exposed_perim=0.0,height=2.5,depth=1.75,window_area=0.0,door_area=0.0,frac_heated=0.0,T_fluid=33,soilk=0.85,watertable_depth=8):
        '''
        Docstring Basesimp

        :int iLoc: Location index 
        :str config: Basesimp configuration identifier. Format BXXX_Y or SXX_Y for basements or slabs, respectively
        :float length: Length of basement. Min of 5 m, max of 30 m
        :float width: Width of basement. Min of 5 m, max of 30 m
        :float overlap: Overlap of insulation. Only used of 'config' has overlap as variable. Default 0.6
        :float rsi_interior: Thermal resistance of interior wall insulation (m2.K/W). Default 0
        :float rsi_exterior: Thermal resistance of exterior wall insulation (m2.K/W). Default 0
        :float rsi_slab: Thermal resistance of slab insulation (m2.K/W)
        :float exposed_perim: Length of perimeter exposed (m). If all exposed, enter 0.
        :float height: Height of basement (m)
        :float depth: Depth of basement (m)
        :float window_area: Total window area (m2). Default 0 m2
        :float door_area: Total door area (m2). Default 0 m2
        :float frac_heated: Fraction of slab heated. If no radiant heating, set to 0. Default 0
        :float T_fluid: Fluid temperature of infloor heating system. 22 to 55 oC
        :float soilk: Soil conductivity (W/m.K). Normal = 0.85, High/moist = 1.275, very wet/permafrost = 1.9
        :float watertable_depth: Water table depth (m). Shallow = 5, Normal = 8, Deep = 12
        '''
        
        sPathDbs = str(pathlib.Path(__file__).resolve().parent).__str__()
        sPathDbs=sPathDbs.replace('\\','/')
        self.__sWthFile = sPathDbs+'/FdnWth.csv'
        sConfigFile = sPathDbs+'/FdnConfigs.csv'
        sConfigFile = sPathDbs+'/FdnConfigs.csv'
        sCFsFile = sPathDbs+'/CornerCFs.csv'
        if not os.path.isfile(self.__sWthFile):
            raise Exception("Cannot find weather database FdnWth.csv. Must be in same folder as class src file.")
        elif not os.path.isfile(sConfigFile):
            raise Exception("Cannot find configurations database FdnConfigs.csv. Must be in same folder as class src file.")
        elif not os.path.isfile(sCFsFile):
            raise Exception("Cannot find corner corrections database CornerCFs.csv. Must be in same folder as class src file.")
        
        self.__soilk:float  = soilk
        self.__wtDepth: float = watertable_depth
        self.__rsiInt: float = rsi_interior
        self.__rsiExt: float = rsi_exterior
        self.__rsiSlab: float = rsi_slab
        self.__loc: int = iLoc
        self.__config: str = config.upper()
        self.__height: float = height
        self.__width: float = width
        self.__length: float = length
        self.__depth: float = depth
        self.__window_area: float=window_area
        self.__door_area: float=door_area
        self.__fracHeated: float = frac_heated
        self.__TRadFluidTemp: float = T_fluid
        self.__exposedPerimLength: float = exposed_perim
        
        self.__checkInputs()
        self.__bIsBasement = True
        if self.__config[0] == 'S':
            self.__setSlabOverride()
            self.__bIsBasement = False
        self.__perim = 2.0*(length+width)
        self.__setConstants()
        self.__setRadiantSlabTemp()
        self.__getWeather(self.__sWthFile) # Get the design weather data
        self.__getConfigs(sConfigFile) # Load the correlation coefficients
        self.__getCornerCF(sCFsFile) # Load the corner correction coefficients
        self.__setOverlap(overlap)
        self.__Wilen = (self.__height-self.__depth)+self.__overlap
        self.__Welen = self.__depth+0.1
        self.__Overlap_06 = self.__overlap/0.6
        self.__setAreas()

        # This method replicates "BS1"
        self.__setBS1()
        
        # This method replicates "Foundation_Calc"
        self.__setFdnCalc()
    
    def getHeatingLoad(self,month=1):
        iDayMonth, iDayYear = self.__getMonthParams(month)
        r1=math.cos(self.__Omega*float(iDayYear)+self.__Ca2)
        if month == 1:
            r2=self.__Ca3
        else:
            iDayMonthP, iDayYearP = self.__getMonthParams(month-1)
            r2=math.cos(self.__Omega*float(iDayYearP)+self.__Ca2)
        
        FHLmon=self.__Ca1*(r2-r1)/iDayMonth
        if self.__bIsBasement:
            TotHeatLoadW=(self.__SagFinal*(self.__TBsmt-self.__Design_H_DBT)*self.__Agfr+FHLmon+self.__alpha2*(self.__TBsmt-self.__DGTemp)+self.__alpha3*(self.__Radiant_Slab_Temperature-self.__DGTemp))*self.__Exposed_Fraction
        else:
            if self.__fracHeated>0.0:
                Tboundary = self.__Radiant_Slab_Temperature
            else:
                Tboundary = self.__TBsmt
            TotHeatLoadW=(self.__SagFinal*(Tboundary-self.__Design_H_DBT)+FHLmon+self.__SbavgFinal*(Tboundary-self.__DGTemp))*self.__Exposed_Fraction
            
        return TotHeatLoadW
    
    def __checkInputs(self):
        if self.__length < 5.0 or self.__length>30.0:
            raise Exception(f"Foundation length must be between 5 and 30 m. Value entered: {self.__length}")
        if self.__width < 5.0 or self.__width>30.0:
            raise Exception(f"Foundation width must be between 5 and 30 m. Value entered: {self.__width}")
        if self.__height < 0.0:
            raise Exception(f"Foundation height must be greater than 0 m. Value entered: {self.__height}")
        if self.__depth < 0.0:
            raise Exception(f"Foundation depth must be greater than 0 m. Value entered: {self.__depth}")
        if self.__depth > self.__height:
            raise Exception(f"Foundation height must be greater than depth. Height of {self.__height} and depth of {self.__depth} entered")
        if self.__TRadFluidTemp < 22.0 or self.__TRadFluidTemp > 50.0:
            raise Exception(f"Hydronic heating fluid temperature must be between 22 and 55. Value enetered is {self.__TRadFluidTemp}")
        # TODO: More error handling

    def __setSlabOverride(self):
        self.__height = 0.0
        self.__depth = 0.05
    
    def __setConstants(self):
        self.__TBsmt = 22.0
        self.__InsValues = {
            'concrete': {
                'walls': 0.116,
                'floors':0.0578
            },
            'wood': {
                'walls': 0.417,
                'floors':0.833
            }
        }
    
    def __setAreas(self):
        self.__Abwag=(self.__height-self.__depth)*2*(self.__length+self.__width)
        self.__Agfr=(self.__Abwag-self.__window_area-self.__door_area)/self.__Abwag
        if self.__exposedPerimLength <= 0:
            self.__Exposed_Fraction = 1.0
        else:
            self.__Exposed_Fraction = min(1,max(self.__exposedPerimLength/self.__perim))

    def __setAlphas(self):
        if int(self.__dfCOeffs['iFndFlag_4']) == 0:
            keyWall = 'concrete'
        else:
            keyWall = 'wood'
        if int(self.__dfCOeffs['iFndFlag_5']) == 0:
            keyFloor = 'concrete'
        else:
            keyFloor = 'wood'
        
        rValWalls = self.__rsiInt+self.__rsiExt+self.__InsValues[keyWall]['walls']
        thickness = rValWalls/13.93
        insideLength = self.__length-2.*thickness
        insidewidth = self.__width-2*thickness
        insideFloorArea  =insideLength*insidewidth
        uWalls = self.__depth*2*(insideLength+insidewidth)/rValWalls
        
        rValFloor = self.__rsiSlab+self.__InsValues[keyFloor]['floors']
        uFloors = insideFloorArea/rValFloor
        
        uVal = max((uWalls+uFloors),0.1)
        if self.__fracHeated > 0.0:
            self.__alpha2 = self.__SbavgFinal*uWalls/uVal
            self.__alpha3 = self.__SbavgFinal*uFloors/uVal
        else:
            self.__alpha2 = self.__SbavgFinal
            self.__alpha3 = 0.0       
    
    def __setRadiantSlabTemp(self):
        if self.__TRadFluidTemp > 0.0:
            self.__Radiant_Slab_Temperature = self.__TBsmt+(self.__fracHeated*(self.__TRadFluidTemp-self.__TBsmt))
        else:
            self.__Radiant_Slab_Temperature = 0.0
        
    def __getWeather(self,sPathDbs):
        lMonths = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
        df = pd.read_csv(sPathDbs, encoding = "ISO-8859-1")
        self.__dfWth = df.loc[df['Seq'] == self.__loc]
        
        self.__DGTemp = float(self.__dfWth['DGTEMP'])
        self.__DegDay =float(self.__dfWth['DegDay'])
        self.__Design_H_DBT = float(self.__dfWth['DHDBT'])
        
        
        self.__AvgAnnTemp=0.0
        Ca1Sum = 0.0
        Ca2Sum = 0.0
        i=1
        for sMon in lMonths:
            self.__AvgAnnTemp+=self.__dfWth[sMon]
            omega =(math.pi/6.0)*(i-0.5)
            i+=1
            Ca1 = self.__dfWth[sMon]*math.sin(omega)
            Ca2 = self.__dfWth[sMon]*math.cos(omega)
            Ca1Sum += Ca1
            Ca2Sum += Ca2
        
        amp0 = -1.0*math.sqrt((Ca1Sum**2.0)+(Ca2Sum**2.0))/6.0
        Ph0 = math.atan(Ca1Sum/Ca2Sum)
        self.__Amps = -1.0*(amp0+0.00197*self.__DegDay-7.8747)
        self.__Phs = Ph0+0.00002128*self.__DegDay-0.0756
    
    def __getConfigs(self,sPathDbs):
        df = pd.read_csv(sPathDbs, encoding = "ISO-8859-1")
        # Find the coefficients for the requested 
        self.__dfCOeffs = df.loc[df['Name'] == self.__config]
        
        # Find the uninsulated configuration coefficients
        iNoIns = int(self.__dfCOeffs['iUnInsul'])
        self.__dfCOeffsNoIns = df.loc[df['#'] == iNoIns]
        
    def __getCornerCF(self,sPathDbs):
        df = pd.read_csv(sPathDbs, encoding = "ISO-8859-1")
        self.__cornerCoeffs = df.values
   
    def __setBS1(self):
        i = 1
        for config in [self.__dfCOeffs,self.__dfCOeffsNoIns]:
            if i == 1:
                rsi1c = self.__getRSI1(config)
                if rsi1c > 1.5:
                    rsi2c = rsi1c
                else:
                    rsi2c = 1.5
            else:
                rsi1c = 0.01
                rsi2c = 0.01           
            
            # Calculate Summo
            rpart1 =(float(config['a1'])+float(config['b1'])*(self.__height-self.__depth)+float(config['cc1'])/self.__soilk)/rsi2c**float(config['d1'])
            rpart2=1/(float(config['e1'])+float(config['i1'])*(self.__overlap**float(config['f1']))*(rsi2c**float(config['g1']))*(self.__height-self.__depth)**float(config['h1']))
            sumuo = rpart1*rpart2+float(config['j1'])
            if i == 1:
                self.__Sag =sumuo*self.__perim
            else:
                self.__SagNoIns =sumuo*self.__perim

            # Calculate Sumur
            rpart1 =(float(config['q2'])+float(config['rr2'])*self.__width)*(float(config['u2'])+float(config['v2'])*self.__soilk)*(float(config['w2'])+float(config['x2'])*self.__depth)
            rpart2 = float(self.__wtDepth**(config['s2']+config['t2']*self.__width+config['y2']*self.__depth))
            rpart3 = float(config['a2']*(self.__depth**config['b2'])*(self.__soilk**config['cc2']))
            rpart4 =float((self.__wtDepth**config['d2'])*(rsi2c**(config['e2']+config['f2']*self.__soilk+config['g2']*self.__depth+config['h2']*self.__overlap)))
            sumur = (rpart1/rpart2)+(rpart3/rpart4)
            
            # Calculate steady-state corner factor
            if self.__depth>2:
                dept = 2.0
            else:
                dept = self.__depth
            if rsi1c>5:
                rss=5.0
            else:
                rss = rsi2c
            if self.__width>10:
                widt=10.0
            else:
                widt = self.__width
            wby2 = widt/2.0
            icol = self.__getIcol(config)
            if icol == 98:
                icol = 3
                rs = 0.0
            else:
                rs = rss
                
            if i == 2: rs = 0.0

            iuse_1=int(2*(icol-1)+1)-1
            Coeffs = self.__cornerCoeffs[iuse_1].tolist()
            r1=Coeffs[0]+Coeffs[1]*rs+Coeffs[2]*self.__soilk+Coeffs[3]*wby2+Coeffs[4]*dept+Coeffs[5]*self.__wtDepth
            r2=Coeffs[6]*rs**2+Coeffs[7]*self.__soilk*rs+Coeffs[8]*wby2*rs+Coeffs[9]*wby2*self.__soilk+Coeffs[10]*wby2**2
            r3=Coeffs[11]*dept*rs+Coeffs[12]*dept*self.__soilk+Coeffs[13]*dept*wby2+Coeffs[14]*dept**2
            r4=Coeffs[15]*self.__wtDepth*rs+Coeffs[16]*self.__wtDepth*self.__soilk+Coeffs[17]*self.__wtDepth*wby2+Coeffs[18]*self.__wtDepth*dept
            Fcs=r1+r2+r3+r4
            thisSbgavg=sumur*(2.0*(self.__length-self.__width)+4*Fcs*self.__width)
            if i == 1:
                self.__Sbgavg=thisSbgavg
            else:
                self.__SbgavgNoIns=thisSbgavg
            
            # Calculate the attenuation
            rpart1=float(config['a3']+config['b3']*self.__soilk+config['cc3']*self.__depth)
            rpart2=float(config['e3']+config['f3']*self.__soilk+config['g3']*self.__depth)
            rpart3=float(rsi2c**(config['h3']+config['i3']*self.__overlap))
            if rpart3>0:
                atten=rpart1+(rpart2/rpart3)
            else:
                atten=rpart1
            
            # Calculate variable corner factor
            iusev=int(2*(icol-1))+2-1
            Coeffs = self.__cornerCoeffs[iusev].tolist()
            r1=float(Coeffs[0]+Coeffs[1]*rs+Coeffs[2]*self.__soilk+Coeffs[3]*wby2+Coeffs[4]*dept+Coeffs[5]*self.__wtDepth)
            r2=float(Coeffs[6]*rs**2+Coeffs[7]*self.__soilk*rs+Coeffs[8]*wby2*rs+Coeffs[9]*wby2*self.__soilk+Coeffs[10]*wby2**2)
            r3=float(Coeffs[11]*dept*rs+Coeffs[12]*dept*self.__soilk+Coeffs[13]*dept*wby2+Coeffs[14]*dept**2)
            r4=float(Coeffs[15]*self.__wtDepth*rs+Coeffs[16]*self.__wtDepth*self.__soilk+Coeffs[17]*self.__wtDepth*wby2+Coeffs[18]*self.__wtDepth*dept)
            Fcv=r1+r2+r3+r4
            
            thisSbgvar=atten*(2*(self.__length-self.__width)+4*self.__width*Fcv)
            if i == 1:
                self.__Sbgvar=thisSbgvar
            else:
                self.__SbgvarNoIns=thisSbgvar
            
            # Calculate Phase
            thisPhase=float(config['a4']+config['b4']/rsi2c**config['cc4'])
            if i == 1:
                self.__fdnPhase=thisPhase
            else:
                self.__fdnPhaseNoIns=thisPhase
            
            i+=1
    
    def __getIcol(self,df):
        # This is a rather convoluted process
        icol1 = int(df['CCF'])
        if icol1 == 99:
            if self.__Overlap_06 <= 0.9999:
                icol1 = 4
            if self.__Overlap_06 > 0.9999 and (self.__Welen/self.__Wilen) > 1.0:
                icol1 = 5
            if self.__Overlap_06 > 0.9999 and (self.__Welen/self.__Wilen) <= 1.0:
                icol1 = 3
        return icol1
    
    def __getRSI1(self,df):
        if self.__bIsBasement:
            iInFlag = self.__getInFlag(df)
            match iInFlag:
                case 1:
                    rsi1=0.0
                case 2:
                    rsi1=max(self.__rsiInt,self.__rsiExt)
                case 3:
                    rsi1=0.88*max(self.__rsiInt,self.__rsiExt)+0.12*self.__rsiSlab
                case 4:
                    rsi1=self.__rsiInt+self.__rsiExt
                case 5:
                    rsi1=(self.__rsiInt+self.__rsiExt)*0.44+0.12*self.__rsiSlab
                case _:
                    raise Exception(f"Insulation flag must be 1 - 5. Value entered is {iInFlag}")
        else:
            if self.__config[2] == 'N':
                return 0
            else:
                return self.__rsiSlab
            
        return rsi1
    
    def __getInFlag(self,df):
        iInFlag = 0
        if int(df['iFndFlag_6']) < 1 and int(df['iFndFlag_7']) < 1:
            iInFlag+=1
        if int(df['iFndFlag_6']) < 1 and int(df['iFndFlag_7']) > 0 and int(df['iFndFlag_8']) < 1:
            iInFlag+=2
        if int(df['iFndFlag_6']) < 1 and int(df['iFndFlag_7']) > 0 and int(df['iFndFlag_8']) > 0:
            iInFlag+=3
        if int(df['iFndFlag_6']) >0 and int(df['iFndFlag_7']) > 0 and int(df['iFndFlag_8']) < 1:
            iInFlag+=4
        if int(df['iFndFlag_6']) >0 and int(df['iFndFlag_7']) <1 and int(df['iFndFlag_8']) < 1:
            iInFlag+=2
        if int(df['iFndFlag_6']) >0 and int(df['iFndFlag_7']) <1 and int(df['iFndFlag_8']) > 0:
            iInFlag+=3
        if int(df['iFndFlag_6']) >0 and int(df['iFndFlag_7']) > 0 and int(df['iFndFlag_8']) > 0:
            iInFlag+=5
        return iInFlag
        
    def __setOverlap(self,entered):
        iId = int(self.__dfCOeffs['#'])
        match iId:
            case 11 | 12 | 116 | 117:
                self.__overlap = entered
            case 93 | 95 | 114 | 115:
                self.__overlap = 0.6
            case 94 | 96:
                self.__overlap = self.__depth-0.2
            case 68 | 69 | 92:
                self.__overlap = self.__depth-0.6
            case _:
                self.__overlap = 0.0

    def __setFdnCalc(self):
        rsi_total=self.__getRSI1(self.__dfCOeffs)
        wrsi= self.__getWRSI()
        denom=math.exp(wrsi*rsi_total)
        if rsi_total<= 0.01:
            self.__SagFinal=self.__SagNoIns
            self.__SbavgFinal=self.__SbgavgNoIns
            self.__SbvarFinal=self.__SbgvarNoIns
            self.__fdnPhaseFinal=self.__fdnPhaseNoIns
        elif rsi_total >= 1.5:
            self.__SagFinal=self.__Sag
            self.__SbavgFinal=self.__Sbgavg
            self.__SbvarFinal=self.__Sbgvar
            self.__fdnPhaseFinal=self.__fdnPhase
        else:
            # Interpolate
            self.__SagFinal=self.__Sag+(self.__SagNoIns-self.__Sag)/denom
            self.__SbavgFinal=self.__Sbgavg+(self.__SbgavgNoIns-self.__Sbgavg)/denom
            self.__SbvarFinal=self.__Sbgvar+(self.__SbgvarNoIns-self.__Sbgvar)/denom
            self.__fdnPhaseFinal=self.__fdnPhase+(self.__fdnPhaseNoIns-self.__fdnPhase)/denom
        
        self.__Omega =2*math.pi/365
        self.__Ca1 =self.__SbvarFinal*self.__Amps/self.__Omega
        self.__Ca2=self.__fdnPhaseFinal-(0.5*math.pi)-self.__Phs
        self.__Ca3=math.cos(self.__Ca2)
        
        # Calculate the "Alphas"
        self.__setAlphas()
    
    def __getWRSI(self):
        if self.__bIsBasement:
            return 2.29
        else:
            return 1.77
    
    def __getMonthParams(self,i):
        match i:
            case 1:
                j=31
                k=31
            case 2:
                j=28
                k=59
            case 3:
                j=31
                k=90
            case 4:
                j=30
                k=120
            case 5:
                j=31
                k=151
            case 6:
                j=30
                k=181
            case 7:
                j=31
                k=212
            case 8:
                j=31
                k=243
            case 9:
                j=30
                k=273
            case 10:
                j=31
                k=304
            case 11:
                j=30
                k=334
            case 12:
                j=31
                k=365
            case _:
                raise Exception(f"Unknown month {i}. Must be 1-12.")
        return j,k