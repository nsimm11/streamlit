import streamlit as st
st.set_page_config(layout="wide")

import pandas as pd
import numpy as np
import math
import os

from pyXSteam.XSteam import XSteam
steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS) # m/kg/sec/°C/bar/W

#inputs
g = 9.81 #m/ss

#PIPING COMPONENTS
konductivity = 16.2 #conductivity of ss
konduct_air = 0.025 #conductivity of air
#inlet
thick_inlet = 0.0095

#header
thick_header = 0.01905 #m

#distribution tubes
thick_tube = 0.0127 #thickness of 3" npt pipe

#condensate drainage line
thick_cond = 0.00762 #thickness of 3/4" npt pipe

#Ambient conditions
pamb = 1 #bar
tamb = 20 #C

#Supplied Steam
psam = 1 #bar

c1, c2, c3 = st.columns(3)
c2.markdown("# Capstone Model")
c2.image("sam.jpg")


st.markdown("### Piping Component Variables")
c1, c2, c3, c4, c5, c6 = st.columns(6)

dia_inlet = float(c1.text_input(label="Inlet Diameter [in]", value=1.75)) * 0.0254

insulated = bool(c2.selectbox(label="Insulated?", options=[False, True]))
thick_insul = float(c2.text_input(label="Insulation Thickness [in]", value=0.5)) * 0.0254

dia_header = float(c3.text_input(label="Header Diameter [in]", value=12)) * 0.0254
len_header = float(c3.text_input(label="Header Length [in] (to dist tube)", value=12)) * 0.0254

dia_tube = float(c4.text_input(label="Distribution Tube Diameter [in]", value=3)) * 0.0254
len_tube = float(c4.text_input(label="Distribution Tube Length [in]", value=10)) * 0.0254

dia_nozzle = float(c5.selectbox(label="Nozzle Diameter [in]", options=[str(1/4), str(1/8), str(1/16)])) * 0.0254
num_nozzles = float(c5.text_input(label="Number of Nozzles", value=str(100)))
dia_cond = float(c6.text_input(label="Condensate Drainage Line Diameter [in]", value=str(0.75))) * 0.0254

st.markdown("### Inlet Steam Conditions")

c1, c2, c3, c4, c5 = st.columns(5)
dryfract_in = float(c1.text_input(label="Dryness Ratio of Inlet Steam", value=0.65))
tsteam = float(c2.text_input(label="Inlet Steam Temperature [C]", value=100))
p_in = float(c3.text_input(label="Inlet Pressure [atm]", value=1))
v_in = float(c4.text_input(label="Inlet Flow Rate [m3/hr]", value=14)) / ((dia_inlet/2)**2 * np.pi * 3600)
z_in = float(c5.text_input(label="Inlet Height [m]", value=0))

st.markdown("### Assumptions")
st.write("""
- Pressure across the SAM-e is atmostpheric, this includes both the nozzles and condensate drainage \n
- The SAM-e is made out of 304SS, with a conductivity constant k of 16.2 W/mk, Radiation and Convection were not considered \n
- The averge height of the nozzles is half the height of the distribution tube, as the nozzles are evenly distributed across them \n
- The height of the inlet and condensate lines are at 0m \n
- The nozzle design completely removes all entrained liquid in the vapour stream so only vapor exits through the nozzles, and only liquid exits through the condensate drainage line
""")

A_in = (dia_header/2)**2 * np.pi

A_nozzle = num_nozzles*(dia_nozzle/2)**2 * np.pi
p_nozzle = pamb
z_nozzle = 0

A_cond = (dia_cond/2)**2 * np.pi
z_cond = 0

#Thermo Calcs
def calc_u(u_x, u_pressure):
  u_h = steamTable.h_px(u_pressure, u_x)
  u = steamTable.u_ph(u_pressure, u_h)
  return u

def calc_sv(sv_x, sv_pressure):
  sv_h = steamTable.h_px(sv_pressure, sv_x)
  sv = steamTable.v_ph(sv_pressure, sv_h)
  return sv

#Fluid Calcs
#mdot from pipe velocity
def calc_mdot(mdot_v_in, mdot_dia, mdot_x, mdot_pressure):
  mdot_sv = calc_sv(mdot_x, mdot_pressure)
  mdot_area = np.pi * (mdot_dia/2)**2
  mdot = (1/mdot_sv) * (mdot_area) * mdot_v_in
  return mdot


#Heat Losses Calcs
def ht_cond(length, t_inner, t_amb, dia_pipe, thickness, konductivity, insulated):
  if insulated:
    cond_losses = 2 * np.pi * length * (t_inner - t_amb) / (((math.log(((dia_pipe/2) + thickness) / (dia_pipe/2))) / konductivity) + ((math.log(((dia_pipe/2) + thickness + thick_insul) / ((dia_pipe/2) + thickness))) / konduct_air))
  else:
    cond_losses = 2 * np.pi * length * (t_inner - t_amb) / ((math.log(((dia_pipe/2) + thickness) / (dia_pipe/2) ))/konductivity)
  return cond_losses

def q_losses(t_in, t_amb, mdot_in, mdot_nozzle, insulated):
  qdot_cond_header = ht_cond(len_header, t_in, t_amb, dia_header, thick_header, konductivity, insulated)*mdot_in #add function here
  qdot_cond_tube = ht_cond(len_tube, t_in, t_amb, dia_tube, thick_tube, konductivity, insulated)*mdot_nozzle #add function here
  qdot_cond = qdot_cond_header + qdot_cond_tube
  qdot_convection = 0 #add function here
  qdot_rad = 0 #add function here
  return qdot_cond + qdot_convection + qdot_rad

#MODEL
def calc_dryfractOut(in_dryfract, pamb, tamb, tsteam, vel_in, dia_in, dia_header, z_in, A_nozzle, A_cond, insulated):
  divi = 10000
  old_error = 10000
  #Assume pressure and temperature are constant throughout SAM-e
  p_in = pamb
  p_nozzle = pamb
  p_cond = pamb #Poor assumption - pipeflow with water and velocity shouldn't be at ambient pressure

  t_in = tsteam
  t_nozzle = tsteam
  t_cond = tsteam
  t_amb = tamb 

  #Inlet Calculations
  mdot_in = calc_mdot(vel_in, dia_in, in_dryfract, p_in)
  u_in = calc_u(in_dryfract, p_in)
  row_in = 1 / (calc_sv(in_dryfract, p_in))
  energy_in = (p_in*101325/row_in)*math.log(p_in*101325)*mdot_in +  mdot_in * (vel_in**2)/2 + g * z_in * mdot_in + u_in*mdot_in
  dryfractguess_out = 1
  error = 1000
  mdot_nozzle = mdot_in * dryfract_in
  mdot_cond = mdot_in * (1 - dryfract_in)
  q_dot = q_losses(t_in, t_amb, mdot_in, mdot_nozzle, insulated)
  count = 0
  while 1:
    count = count + 1
    #Outlet Calculations (Conditions at Nozzle assumed to be 100% dry steam, condensate outlet assumed to be 100% sat liquid)
    p_out = p_in
    uguess_nozzle = calc_u(dryfractguess_out, p_nozzle)
    uguess_cond = calc_u((0), p_nozzle)
    row_nozzle = 1 / calc_sv((dryfractguess_out), p_nozzle)
    row_cond =  1 / calc_sv(0, p_nozzle)  #100% liquid
    vel_nozzle = mdot_nozzle / (A_nozzle * row_nozzle)
    vel_cond = mdot_cond / (A_cond * row_cond)
    energy_nozzle = (p_nozzle*101325/row_nozzle)*math.log(p_nozzle*101325)*mdot_nozzle +  mdot_nozzle*(vel_nozzle**2)/2 + g*z_nozzle*mdot_nozzle + uguess_nozzle*mdot_nozzle
    energy_cond = (p_cond*101325/row_cond)*math.log(p_cond*101325)*mdot_cond +  mdot_cond*(vel_cond**2)/2 + g*z_cond*mdot_cond + uguess_cond*mdot_cond
    error = energy_in - (energy_nozzle + energy_cond) - q_dot
    print(error, dryfractguess_out)
    if abs(error) < 0.1:
      break
    if abs(old_error) < abs(error):
        print("ERROR GOING UP")
        divi = divi * 10

    dryfractguess_out = dryfractguess_out + (error / divi)
    if dryfractguess_out > 1:
      st.markdown("# DRY FRACTION GREATER THAN 1, ERROR:" + str(error))
      break
      
    old_error = error
    if count > 1000:
        st.markdown("# COULD NOT CONVERGE, ERROR: " + str(error))
        break
    
  dryfract_out = dryfractguess_out
  return dryfract_out, q_dot

button = st.button("Calculate")

if button:
    dfo, qd = calc_dryfractOut(dryfract_in, pamb, tamb, tsteam, v_in, dia_inlet, dia_header, z_in, A_nozzle, A_cond, insulated)
    st.session_state["qdot"] = qd
    st.session_state["liqfract"] = dfo

    st.markdown("### Results")
    st.write("Heat Losses per second: " + str(np.round(st.session_state['qdot'],2)) + " KJ/s")
    st.write("Vapor Percentage of Outgoing Steam: " + str(np.round(st.session_state['liqfract'],4) * 100) + " %")

