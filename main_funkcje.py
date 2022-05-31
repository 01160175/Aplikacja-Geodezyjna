# -*- coding: utf-8 -*-
"""
Created on Sun May 22 12:41:34 2022

@author: zmigr
"""
import math as m
import numpy as np
import sys
from PyQt5.QtWidgets import QDialog, QApplication, QMessageBox
import PyQt5.QtGui
from App_geod import*

class MyForm(QDialog):
    def __init__(self): 
        super().__init__()
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        self.setWindowIcon(QtGui.QIcon('ikona.png'))
        self.setWindowTitle("Aplikacja Geodezyjna")
        
        self.ui.p2000.clicked.connect(self.przelicz2000)
        self.ui.p1992.clicked.connect(self.przelicz1992)
        self.ui.XYZ.clicked.connect(self.przeliczXYZ)
        
        self.ui.p2000_st.clicked.connect(self.przelicz2000_tab2)
        self.ui.p1992_st.clicked.connect(self.przelicz1992_tab2)
        self.ui.XYZ_2.clicked.connect(self.przeliczXYZ_tab2)

        
        self.ui.WGS84.clicked.connect(self.elipsoida)
        self.ui.GRS80.clicked.connect(self.elipsoida)
        
        self.ui.WGS84_st.clicked.connect(self.elipsoida)
        self.ui.GRS80_st.clicked.connect(self.elipsoida)
        
        self.ui.str5.toggled.connect(self.strefy)
        self.ui.str6.toggled.connect(self.strefy)
        self.ui.str7.toggled.connect(self.strefy)
        self.ui.str8.toggled.connect(self.strefy)
        
        self.ui.str5_st.toggled.connect(self.strefy2)
        self.ui.str6_st.toggled.connect(self.strefy2)
        self.ui.str7_st.toggled.connect(self.strefy2)
        self.ui.str8_st.toggled.connect(self.strefy2)
        
        
        self.show()
        
        
# STREFY ---------------------------------------------------------------------------------------

        """   
        
        Funkcje dobierające parametry strefy dla układu 2000.
        
        """
        
    def strefy(self):
        if self.ui.str5.toggled:
            self.s = 5
            self.lam0 = 15
            self.lam0 = np.radians(self.lam0)
                
        elif self.ui.str6.toggled:
            self.s = 6
            self.lam0 = 18
            self.lam0 = np.radians(self.lam0)
                
        elif self.ui.str7.toggled:
            self.s = 7
            self.lam0 = 21
            self.lam0 = np.radians(self.lam0)
                
        elif self.ui.str8.toggled:
            self.s = 8
            self.lam0 = 24
            self.lam0 = np.radians(self.lam0)
           
        return self.s,self.lam0
    
    def strefy2(self):
        if self.ui.str5_st.toggled:
            self.s = 5
            self.lam0 = 15
            self.lam0 = np.radians(self.lam0)
                
        elif self.ui.str6_st.toggled:
            self.s = 6
            self.lam0 = 18
            self.lam0 = np.radians(self.lam0)
                
        elif self.ui.str7_st.toggled:
            self.s = 7
            self.lam0 = 21
            self.lam0 = np.radians(self.lam0)
                
        elif self.ui.str8_st.toggled:
            self.s = 8
            self.lam0 = 24
            self.lam0 = np.radians(self.lam0)
           
        return self.s,self.lam0
    
# Uklady 2000 --------------------------------------------------------------------------------------- 
    
        """   
        Funkcje przeliczające współrzędne geodezyjne na współrzędne układu 2000.
        
        Parametry
        -------
        fi  [float] : szerokość geodezyjna [rad]
        lam [float] : długość geodezyjna [rad]
        a   [float] : dłuższa półoś elipsoidy [m]
        e2  [float] : mimośród elipsoidy
        
        Wyniki
        -------
        x00 [float] : współrzędna w układzie 2000 [m]
        y00 [float] : współrzędna w układzie 2000 [m]
        
        """

    def uk00(self,fi,lam):
        m = 0.999923
    
        N = self.a/np.sqrt(1-self.ecc*np.sin(fi)**2)
        e2p = self.ecc/(1-self.ecc)
        t = np.arctan(fi)
        n2 = e2p * (np.cos(fi))**2
        lam = np.degrees(lam)
    
        lam = np.radians(lam)
        
        l = lam - self.lam0
    
        A0 = 1 - (self.ecc/4) - ((3*(self.ecc**2))/64) - ((5*(self.ecc**3))/256)
        A2 = (3/8) * (self.ecc + ((self.ecc**2)/4) + ((15 * (self.ecc**3))/128))
        A4 = (15/256) * (self.ecc**2 + ((3*(self.ecc**3))/4))
        A6 = (35 * (self.ecc**3))/3072
    
        sig = self.a * ((A0*fi) - (A2*np.sin(2*fi)) +
                   (A4*np.sin(4*fi)) - (A6*np.sin(6*fi)))
        x = sig + ((l**2)/2) * (N*np.sin(fi)*np.cos(fi)) * (1 + ((l**2)/12) * ((np.cos(fi))**2) * (5 - t**2 + 9*n2 +
                                                                                                   4*(n2**2)) + ((l**4)/360) * ((np.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 * (t**2))))
        y = l * (N*np.cos(fi)) * (1 + ((((l**2)/6) * (np.cos(fi))**2) * (1-(t**2) + n2)) +
                                  (((l**4)/(120)) * (np.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
    
        x00 = m * x
        y00 = m * y + (self.s*1000000) + 500000
        
        self.ui.wynik10.setText('X = ' + str(round(x00,3)) + ', ' + 'Y = ' + str(round(y00,3)))
        
        
        
    def uk00_tab2(self,fi,lam):
        m = 0.999923
    
        N = self.a/np.sqrt(1-self.ecc*np.sin(fi)**2)
        e2p = self.ecc/(1-self.ecc)
        t = np.arctan(fi)
        n2 = e2p * (np.cos(fi))**2
        lam = np.degrees(lam)
    
        lam = np.radians(lam)
        
        l = lam - self.lam0
    
        A0 = 1 - (self.ecc/4) - ((3*(self.ecc**2))/64) - ((5*(self.ecc**3))/256)
        A2 = (3/8) * (self.ecc + ((self.ecc**2)/4) + ((15 * (self.ecc**3))/128))
        A4 = (15/256) * (self.ecc**2 + ((3*(self.ecc**3))/4))
        A6 = (35 * (self.ecc**3))/3072
    
        sig = self.a * ((A0*fi) - (A2*np.sin(2*fi)) +
                   (A4*np.sin(4*fi)) - (A6*np.sin(6*fi)))
        x = sig + ((l**2)/2) * (N*np.sin(fi)*np.cos(fi)) * (1 + ((l**2)/12) * ((np.cos(fi))**2) * (5 - t**2 + 9*n2 +
                                                                                                   4*(n2**2)) + ((l**4)/360) * ((np.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 * (t**2))))
        y = l * (N*np.cos(fi)) * (1 + ((((l**2)/6) * (np.cos(fi))**2) * (1-(t**2) + n2)) +
                                  (((l**4)/(120)) * (np.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
    
        x00 = m * x
        y00 = m * y + (self.s*1000000) + 500000
        
        self.ui.wynikst.setText('X = ' + str(round(x00,3)) + ', ' + 'Y = ' + str(round(y00,3)))
        
        
# Uklady 1992 ---------------------------------------------------------------------------------------   
        """   
        Funkcje przeliczające współrzędne geodezyjne na współrzędne układu 1992.
        
        Parametry
        -------
        fi  [float] : szerokość geodezyjna [rad]
        lam [float] : długość geodezyjna [rad]
        a   [float] : dłuższa półoś elipsoidy [m]
        e2  [float] : mimośród elipsoidy
        
        Wyniki
        -------
        x92 [float] : współrzędna w układzie 1992 [m]
        y92 [float] : współrzędna w układzie 1992 [m]  
        
        """
        
    def uk92(self):
        lam0 = np.radians(19)
        b2 = (self.a**2)*(1 - self.ecc)
        ep2 = ((self.a**2 - b2))/(b2)
        t = np.tan(self.fi)
        n2 = ep2 * ((np.cos(self.fi))**2)
        N =  self.a / (np.sqrt(1 - self.ecc * (np.sin(self.fi)) ** 2))
        A0 = 1-(self.ecc/4) - ((3*(self.ecc**2))/64)-((5*(self.ecc**3))/256)
        A2 = (3/8)*(self.ecc+((self.ecc**2)/4)+((15*(self.ecc**3))/128))
        A4 = (15/256)*((self.ecc**2)+((3*(self.ecc**3))/4))
        A6 = (35*(self.ecc**3))/3072
        dlam = self.lam - lam0
        sigma = self.a * (A0 * (self.fi) - (A2 * np.sin(2 * self.fi)) + (A4 * np.sin(4 * self.fi)) - (A6 * np.sin(6 * self.fi)))
        x = sigma + ((dlam**2)/2) * N * np.sin(self.fi) * np.cos(self.fi) * (1 + ((dlam**2)/12) * ((np.cos(self.fi))**2) * (5 - t**2 + 9 * n2 + 4 * (n2**2)) + ((dlam**4)/360) * ((np.cos(self.fi))**4) * (61 - (58 * (t**2)) + (t**4) + (270 * n2) - (330 * n2 * (t**2))))
        y = dlam * (N * np.cos(self.fi)) * (1 + ((((dlam**2)/6) * (np.cos(self.fi))**2) * (1 - t**2 + n2)) + (((dlam**4)/(120)) * (np.cos(self.fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14 * n2) - (58 * n2 * (t**2))))
        m = 0.9993
        x92 = x * m - 5300000
        y92 = y * m + 500000
    
        self.ui.wynik92.setText('X = ' + str(round(x92,3)) + ', ' + 'Y = ' + str(round(y92,3)))
        
    def uk92_tab2(self):
        lam0 = np.radians(19)
        b2 = (self.a**2)*(1 - self.ecc)
        ep2 = ((self.a**2 - b2))/(b2)
        t = np.tan(self.fi)
        n2 = ep2 * ((np.cos(self.fi))**2)
        N =  self.a / (np.sqrt(1 - self.ecc * (np.sin(self.fi)) ** 2))
        A0 = 1-(self.ecc/4) - ((3*(self.ecc**2))/64)-((5*(self.ecc**3))/256)
        A2 = (3/8)*(self.ecc+((self.ecc**2)/4)+((15*(self.ecc**3))/128))
        A4 = (15/256)*((self.ecc**2)+((3*(self.ecc**3))/4))
        A6 = (35*(self.ecc**3))/3072
        dlam = self.lam - lam0
        sigma = self.a * (A0 * (self.fi) - (A2 * np.sin(2 * self.fi)) + (A4 * np.sin(4 * self.fi)) - (A6 * np.sin(6 * self.fi)))
        x = sigma + ((dlam**2)/2) * N * np.sin(self.fi) * np.cos(self.fi) * (1 + ((dlam**2)/12) * ((np.cos(self.fi))**2) * (5 - t**2 + 9 * n2 + 4 * (n2**2)) + ((dlam**4)/360) * ((np.cos(self.fi))**4) * (61 - (58 * (t**2)) + (t**4) + (270 * n2) - (330 * n2 * (t**2))))
        y = dlam * (N * np.cos(self.fi)) * (1 + ((((dlam**2)/6) * (np.cos(self.fi))**2) * (1 - t**2 + n2)) + (((dlam**4)/(120)) * (np.cos(self.fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14 * n2) - (58 * n2 * (t**2))))
        m = 0.9993
        x92 = x * m - 5300000
        y92 = y * m + 500000
    
        self.ui.wynik92_2.setText('X = ' + str(round(x92,3)) + ', ' + 'Y = ' + str(round(y92,3)))
        
     
# flh - XYZ ---------------------------------------------------------------------------------------      
     
        """
        Funkcje przeliczające współrzędne geodezyjne  
        na współrzędne geocentryczne (ECEF).

        Parametry
        ----------
        fi  [float] : szerokość geodezyjna [rad]
        lam [float] : długość geodezyjna [rad]
        h   [float] : wysokość elipsoidalna [m]
        a   [float] : dłuższa półoś elipsoidy [m]
        e2  [float] : mimośród elipsoidy
        
        Wyniki
        -------
        X   [float] : współrzędna geocentryczna (ECEF) [m]
        Y   [float] : współrzędna geocentryczna (ECEF) [m]
        Z   [float] : współrzędna geocentryczna (ECEF) [m]
    
        """

    def filh_2_XYZ(self):
        N = self.a/m.sqrt(1-self.ecc*m.sin(self.fi)**2)
        X = (N + self.h) * m.cos(self.fi) * m.cos(self.lam)
        Y = (N + self.h) * m.cos(self.fi) * m.sin(self.lam)
        Z = (N*(1-self.ecc) + self.h) * m.sin(self.fi)
        
        self.ui.wynikX.setText('X = ' + str(X))
        self.ui.wynikY.setText('Y = ' + str(Y))
        self.ui.wynikZ.setText('Z = ' + str(Z))
        
    def filh_2_XYZ_tab2(self):
        N = self.a/m.sqrt(1-self.ecc*m.sin(self.fi)**2)
        X = (N + self.h) * m.cos(self.fi) * m.cos(self.lam)
        Y = (N + self.h) * m.cos(self.fi) * m.sin(self.lam)
        Z = (N*(1-self.ecc) + self.h) * m.sin(self.fi)
        
        self.ui.wynikX_2.setText('X = ' + str(X))
        self.ui.wynikY_2.setText('Y = ' + str(Y))
        self.ui.wynikZ_2.setText('Z = ' + str(Z))
        
# ELIPSOIDA ---------------------------------------------------------------------------------------
    
    """
    
    Funkcja definiująca parametry elipsoidy w zależnosci od zaznaczonego "radioButtona"
    
    """
    def elipsoida(self):
        if self.ui.WGS84.isChecked()==True:
            self.a = 6378137.0
            self.b = 6356752.31424518
            self.flattening = (self.a-self.b)/ self.a
            self.ecc = 2*self.flattening - self.flattening**2
            
        else:
            self.a = 6378137.0
            self.b = 6356752.31414036
            self.flattening = (self.a-self.b)/ self.a
            self.ecc = 2*self.flattening - self.flattening**2

        
# TAB 1 --------------------------------------------------------------------------------------
    
    """
    
    Funkcje wykonujące przeliczenia dla pierwszej zakładki (Stopnie dziesiętne)
    
    fi i lam przeliczane są do radianów
    
    """
    
    def przelicz2000(self):
        if len(self.ui.FI.text()) != 0 and len(self.ui.LAM.text()) != 0:
            fi = (float(self.ui.FI.text())) * (m.pi/180)
            lam = (float(self.ui.LAM.text())) * (m.pi/180)
            
        self.uk00(fi,lam)
        

    def przelicz1992(self):
        if len(self.ui.FI.text()) != 0 and len(self.ui.LAM.text()) != 0:
            self.fi = (float(self.ui.FI.text())) * (m.pi/180)
            self.lam = (float(self.ui.LAM.text())) * (m.pi/180)
            
        self.uk92()

        
    def przeliczXYZ(self):
        if len(self.ui.FI.text()) != 0 and len(self.ui.LAM.text()) != 0 and len(self.ui.Htab1.text()) != 0:
            self.fi = (float(self.ui.FI.text())) * (m.pi/180)
            self.lam = (float(self.ui.LAM.text())) * (m.pi/180)
            self.h = float(self.ui.Htab1.text())
            
        self.filh_2_XYZ()
    
            
# TAB 2 --------------------------------------------------------------------------------------
   
    """
    
    Funkcje wykonujące przeliczenia dla drugiej zakładki (Stopnie)
    
    fi i lam przeliczane są do radianów najpierw będąc konwertowane do stopni dziesiętnych
    
    """
     
    def przelicz2000_tab2(self):
        if len(self.ui.fist.text()) != 0 and len(self.ui.fimin.text()) != 0 and len(self.ui.fisek.text()) != 0 and len(self.ui.lamst.text()) != 0 and len(self.ui.lammin.text()) != 0 and len(self.ui.lamsek.text()) != 0:
            fi = (float(self.ui.fist.text()) + (float(self.ui.fimin.text()))/60 + (float(self.ui.fisek.text()))/3600) * (m.pi/180)
            lam = (float(self.ui.lamst.text()) + (float(self.ui.lammin.text()))/60 + (float(self.ui.lamsek.text()))/3600) * (m.pi/180)
            
        self.uk00_tab2(fi,lam)
        
    
    
    def przelicz1992_tab2(self):
        if len(self.ui.fist.text()) != 0 and len(self.ui.fimin.text()) != 0 and len(self.ui.fisek.text()) != 0 and len(self.ui.lamst.text()) != 0 and len(self.ui.lammin.text()) != 0 and len(self.ui.lamsek.text()) != 0:
            self.fi = (float(self.ui.fist.text()) + (float(self.ui.fimin.text()))/60 + (float(self.ui.fisek.text()))/3600) * (m.pi/180)
            self.lam = (float(self.ui.lamst.text()) + (float(self.ui.lammin.text()))/60 + (float(self.ui.lamsek.text()))/3600) * (m.pi/180)
            
        self.uk92_tab2()
        
        
    def przeliczXYZ_tab2(self):
        if len(self.ui.fist.text()) != 0 and len(self.ui.fimin.text()) != 0 and len(self.ui.fisek.text()) != 0 and len(self.ui.lamst.text()) != 0 and len(self.ui.lammin.text()) != 0 and len(self.ui.lamsek.text()) != 0 and len(self.ui.Htab2.text()) != 0:
            self.fi = (float(self.ui.fist.text()) + (float(self.ui.fimin.text()))/60 + (float(self.ui.fisek.text()))/3600) * (m.pi/180)
            self.lam = (float(self.ui.lamst.text()) + (float(self.ui.lammin.text()))/60 + (float(self.ui.lamsek.text()))/3600) * (m.pi/180)
            self.h = float(self.ui.Htab2.text())
            
        self.filh_2_XYZ_tab2()
        
if __name__=="__main__":
    app = QApplication(sys.argv)
    w = MyForm()
    w.show()
    sys.exit(app.exec_())