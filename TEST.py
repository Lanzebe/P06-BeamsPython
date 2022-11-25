from BeamLib import *

BeamStructure = Beam(1)

Shigley = BeamLab(BeamStructure)
Shigley.AddLoadCase('Shigley_3')
Shigley.LoadCases['Shigley_3'].AddFixedSupport('a', Position = 0)
#Shigley.LoadCases['Shigley_3'].AddDistributedLoading('b', Positions = [0,0.5],Magnitudes = [-1,-1])
#Shigley.LoadCases['Shigley_3'].AddDistributedLoading('c', Positions = [0.4,1],Magnitudes = [-1,-1])

#Shigley.LoadCases['Shigley_3'].AddDistributedLoading('c', Positions = [0.1,0.2],Magnitudes = [-4,-1])
#Shigley.LoadCases['Shigley_3'].AddDistributedLoading('d', Positions = [0.3,0.5],Magnitudes = [-4,-1])
Shigley.LoadCases['Shigley_3'].AddDistributedLoading('e', Positions = [0.1,0.9],Magnitudes = [-4,-4])

#Shigley.LoadCases['Shigley_3'].AddForce('b', Position = 0.5 ,Magnitude = -1)
#Shigley.LoadCases['Shigley_3'].AddForce('c', Position = 0.6 ,Magnitude = -0.2)


Shigley.SolveAll()
Shigley.PlotAll()
