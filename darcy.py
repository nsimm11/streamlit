import numpy as np

def calc_df(flow, diameter, roughness, density, viscosity):

    """
    flow - m3/hr
    diameter - m
    roughness - mm
    density - kg/m3
    visc - Pa.s
    """

    #Turn all values into floating points
    flow = float(flow)/3.6 #into l/s
    diameter = float(diameter)
    roughness = float(roughness)
    density = float(density)
    viscosity = float(viscosity)

    pipearea = float(diameter**2 / 4 * np.pi) #Pipe Area in m^2
    velocity = (flow/1000)/pipearea #Flow Velocity in m/s
    reynolds = velocity*diameter*density/viscosity # Reynolds Number for Full Circular Pipe

    friction = 0.08 #Starting Friction Factor
    while 1:
        leftF = 1 / friction**0.5 #Solve Left side of Eqn
        rightF = - 2 * np.log10(2.51/(reynolds * friction**0.5)+(roughness/1000)/(3.72*diameter)) # Solve Right side of Eqn
        friction = friction - 0.000001 #Change Friction Factor
      #  print(leftF)
      #  print(rightF)
      #  print(friction)
        if (rightF - leftF <= 0): #Check if Left = Right
            break
    return friction
