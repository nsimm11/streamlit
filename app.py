import re
import streamlit as st
st.set_page_config(layout="wide")

import pandas as pd
import numpy as np

from pyXSteam.XSteam import XSteam
steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS) # m/kg/sec/°C/bar/W

from darcy import calc_df

#inputs
g = 9.81 #m/ss

#PIPING COMPONENTS
konductivity = 16.2 #conductivity of ss
konduct_air = 0.025 #conductivity of air

roughness = 0.015

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

opts_headerDia = [12,20]
opts_inletDia = [1, 1.5, 2, 3]
opts_tubeDia = [1, 1.5, 2, 3]
opts_condDia = [0.5, 0.75, 1, 1.5, 2, 3]



st.markdown("### Piping Component Variables")
c1, c2, c3, c4, c5, c6 = st.columns(6)

dia_inlet = float(c1.selectbox(label="Inlet Diameter [in]", options=opts_inletDia)) * 0.0254

insulated = bool(c2.selectbox(label="Insulated?", options=[False, True]))
thick_insul = float(c2.text_input(label="Insulation Thickness [in]", value=0.25)) * 0.0254

dia_header = float(c3.selectbox(label="Header Diameter [in]", options=opts_headerDia)) * 0.0254
len_header = float(c3.text_input(label="Header Length [in] (to dist tube)", value=12)) * 0.0254

dia_tube = float(c4.selectbox(label="Distribution Tube Diameter [in]", options=opts_tubeDia)) * 0.0254
len_tube = float(c4.text_input(label="Distribution Tube Length [in]", value=10)) * 0.0254

dia_nozzle = float(c5.selectbox(label="Nozzle Diameter [in]", options=[str(1/4), str(1/8), str(1/16)])) * 0.0254
num_nozzles = float(c5.text_input(label="Number of Nozzles", value=str(100)))

dia_cond = float(c6.selectbox(label="Condensate Drainage Line Diameter [in]", options=opts_condDia)) * 0.0254

st.markdown("### Inlet Steam Conditions")

c1, c2, c3, c4 = st.columns(4)
dryfract_in = float(c1.text_input(label="Dryness Ratio of Inlet Steam", value=0.65))
tsteam = float(c2.text_input(label="Inlet Steam Temperature [C]", value=100))
p_in = float(c3.text_input(label="Inlet Pressure [atm]", value=1))
v_in = float(c4.text_input(label="Inlet Flow Rate [m3/hr]", value=14)) / ((dia_inlet/2)**2 * np.pi * 3600)
z_in = 0

st.markdown("### Assumptions")
st.write("""
- Pressure across the SAM-e is atmostpheric, this includes both the nozzles and condensate drainage \n
- The SAM-e is made out of 304SS, with a conductivity constant k of 16.2 W/mk, Radiation and Convection were not considered \n
- Minor Line Losses were not considered \n
- The averge height of the nozzles is half the height of the distribution tube, as the nozzles are evenly distributed across them \n
- The height of the inlet and condensate lines are at 0m \n
- The nozzle design completely removes all entrained liquid in the vapour stream so only vapor exits through the nozzles, and only liquid exits through the condensate drainage line
- Does not consider superheated vapor or supercooled liquid
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

def calc_visc(visc_x, visc_pressure):
  visc_h = steamTable.h_px(visc_pressure, visc_x)
  visc = steamTable.my_pt(visc_pressure, 100)
  return visc

#Fluid Calcs
#mdot from pipe velocity
def calc_mdot(mdot_v_in, mdot_dia, mdot_x, mdot_pressure):
  mdot_sv = calc_sv(mdot_x, mdot_pressure)
  mdot_area = np.pi * (mdot_dia/2)**2
  mdot = (1/mdot_sv) * (mdot_area) * mdot_v_in
  return mdot

def calc_reynolds():
  re = 0
  return re

def calc_darcy(flow, diameter, density, viscosity):
  darcy = calc_df(flow, diameter, roughness, density, viscosity)
  return darcy

def calc_majorLosses(pipes, density, viscosity):
  majL = 0
  for pipe in pipes:
    length = pipe[0]
    vel = pipe[1]
    dia = pipe[2]
    darcy = calc_darcy(vel * np.pi*(dia/2)**2, dia, density, viscosity)
    majL_pipe = (darcy * length * vel**2) / (dia * 2 * g)
    majL = majL + majL_pipe
    print("MAJOR: ", majL)
  return majL

def calc_lossCoefficent():

  K = 0

  return K

def calc_minorLosses():

  minL = 0

  return minL

def calc_pipeLosses(f1, f2, l1, l2, d1, d2, density, viscosity):
  #Losses across various pipes
  pipe1 = [l1, f1/(np.pi*(d1/2)**2), d1] # length, vel, dia
  pipe2 = [l2, f2/(np.pi*(d2/2)**2), d2] # length, vel, dia

  mjL = calc_majorLosses([pipe1, pipe2], density, viscosity)

  mnL = calc_minorLosses()
  return mjL + mnL



#Heat Losses Calcs
def ht_cond(length, t_inner, t_amb, dia_pipe, thickness, konductivity, insulated):
  if insulated:
    cond_losses = 2 * np.pi * length * (t_inner - t_amb) / (((np.log(((dia_pipe/2) + thickness) / (dia_pipe/2))) / konductivity) + ((np.log(((dia_pipe/2) + thickness + thick_insul) / ((dia_pipe/2) + thickness))) / konduct_air))
  else:
    cond_losses = 2 * np.pi * length * (t_inner - t_amb) / ((np.log(((dia_pipe/2) + thickness) / (dia_pipe/2) ))/konductivity)
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
  visc_in = calc_visc(in_dryfract, p_in)
  energy_in = (p_in*101325/row_in)*np.log(p_in*101325)*mdot_in +  mdot_in * (vel_in**2)/2 + g * z_in * mdot_in + u_in*mdot_in
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

    pipeLosses = calc_pipeLosses((mdot_in / row_in), (mdot_nozzle / row_nozzle), len_header, len_tube, dia_header, dia_tube, (1/row_in), visc_in) * mdot_in * g

    energy_nozzle = (p_nozzle*101325/row_nozzle)*np.log(p_nozzle*101325)*mdot_nozzle +  mdot_nozzle*(vel_nozzle**2)/2 + g*z_nozzle*mdot_nozzle + uguess_nozzle*mdot_nozzle
    energy_cond = (p_cond*101325/row_cond)*np.log(p_cond*101325)*mdot_cond +  mdot_cond*(vel_cond**2)/2 + g*z_cond*mdot_cond + uguess_cond*mdot_cond
    error = energy_in - (energy_nozzle + energy_cond) - q_dot - pipeLosses
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
  return dryfract_out, q_dot, vel_nozzle, pipeLosses

button = st.button("Calculate")

if button:
    dfo, qd, v_out, pipelosses = calc_dryfractOut(dryfract_in, pamb, tamb, tsteam, v_in, dia_inlet, dia_header, z_in, A_nozzle, A_cond, insulated)
    st.session_state["qdot"] = qd
    st.session_state["liqfract"] = dfo

    st.markdown("## Results")
    st.subheader("Heat Losses per second: " + str(np.round(st.session_state['qdot'],2)) + " KJ/s")
    st.subheader("Pipe Losses per second: " + str(round(pipelosses,4)) + (" KJ/s"))
    st.subheader("Vapor Percentage of Outgoing Steam: " + str(np.round(st.session_state['liqfract'],4) * 100) + " %")
    st.subheader("Nozzle Steam Exit Velocity: " + str(np.round(v_out,4)) + " m/s")


