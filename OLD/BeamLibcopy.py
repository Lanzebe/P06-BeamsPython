import sympy as sp
from math import radians
from sympy.plotting import plot

class BeamLab():
    def __init__(self,Beam) -> None:
        self.Beam = Beam
        self.LoadCases = []
        self.xdim = sp.var('x')

    def AddLoadCase(self,LoadCaseName):
        self.LoadCases.append(LoadCase(LoadCaseName,self))

class Beam():
    def __init__(self,Lenght,ElasticModulus = 210E9,I_Value= (4E6)*1E-12) -> None:
        self.Lenght = Lenght
        self.ElasticModulus = ElasticModulus    #[Pa]
        self.I_Value = I_Value                  #[m^4]

        print(ElasticModulus)
        print(I_Value)

class LoadCase():
    def __init__(self,LoadCaseName,LabEnvironment) -> None:
        self.LoadCaseName = LoadCaseName
        self.Forces = []
        self.Moments = []
        self.FixedSupports = []
        self.PinnedSupports = []
        self.VerticalSliderSupports = []
        self.Lab = LabEnvironment

        self.BoundaryConditions = [ShearForceBoundaryCondition(self.Lab.Beam.Lenght),BendingMomentBoundaryCondition(self.Lab.Beam.Lenght)]
        self.EQ1 = 0
        self.EQ2 = 0
        self.EQ3 = 0
        self.EQ4 = 0
        self.EQ5 = 0
        self.SolvableSetOfEquations = []
        self.C1 = sp.symbols('C_1')
        self.C2 = sp.symbols('C_2')
        self.UnknownVariables = [self.C1,self.C2]
    ##### Building of Loadcase:
    def AddForce(self,ForceName,Magnitude,Position):
        self.Forces.append(Force(ForceName,Magnitude,Position,self))

    def AddMoment(self,MomentName,Magnitude,Position):
        self.Moments.append(Moment(MomentName,Magnitude,Position,self))

    def AddFixedSupport(self,FixedSupportName,Position):
        self.FixedSupports.append(FixedSupport(FixedSupportName,Position,self))
    
    def AddPinnedSupport(self,PinnedSupportName,Position):
        self.PinnedSupports.append(PinnedSupport(PinnedSupportName,Position,self))
    
    def AddVerticalSliderSupport(self,VerticalSliderSupportName,Position):
        self.VerticalSliderSupports.append(VerticalSliderSupport(VerticalSliderSupportName,Position,self))
    ##### Building of Equations:
    def BuildEQ1(self):
        self.EQ1 = 0

    def BuildEQ2(self):
        #self.EQ2 = sp.integrate(self.EQ1,(self.Lab.xdim,0,self.Lab.Beam.Lenght))
        self.EQ2 = sp.integrate(self.EQ1,self.Lab.xdim)
        for i in self.Forces:
            self.EQ2 += i.Equation

    def BuildEQ3(self):
        self.EQ3 = sp.integrate(self.EQ2,self.Lab.xdim)
        for i in self.Moments:
            self.EQ3 += i.Equation

    def BuildEQ4(self):
        self.EQ4 = sp.integrate(self.EQ3/(self.Lab.Beam.ElasticModulus*self.Lab.Beam.I_Value),self.Lab.xdim) + self.C1
        print(self.EQ4)
        
    def BuildEQ5(self):
        self.EQ5 = sp.integrate(self.EQ4,self.Lab.xdim) + self.C2

    def BuildEquations(self):
        self.BuildEQ1()
        self.BuildEQ2()
        self.BuildEQ3()
        self.BuildEQ4()
        self.BuildEQ5()

    def ExtractSolvableSetOFEquations(self):
        for i in self.BoundaryConditions:
            if type(i) == ShearForceBoundaryCondition:
                self.SolvableSetOfEquations.append(self.EQ2.subs({'x':i.Position})-i.Value)
            if type(i) == BendingMomentBoundaryCondition:
                self.SolvableSetOfEquations.append(self.EQ3.subs({'x':i.Position})-i.Value)
            if type(i) == AngleBoundaryCondition:
                self.SolvableSetOfEquations.append(self.EQ4.subs({'x':i.Position})-i.Value)
            if type(i) == DisplacementBoundaryCondition:
                self.SolvableSetOfEquations.append(self.EQ5.subs({'x':i.Position})-i.Value)

    def SolveLoadCase(self):
        self.BuildEquations()
        self.ExtractSolvableSetOFEquations()

        Sol = list(sp.linsolve(self.SolvableSetOfEquations,self.UnknownVariables))[0]

        dictionary = dict(zip(self.UnknownVariables,Sol))
        print("this is the solved values",dictionary)
        ShearForce = self.EQ2.subs(dictionary)
        BendingMoment = self.EQ3.subs(dictionary)
        Angle = self.EQ4.subs(dictionary)
        Displacement = self.EQ5.subs(dictionary)
        
        plot(ShearForce,(self.Lab.xdim,0,self.Lab.Beam.Lenght))
        plot(BendingMoment,(self.Lab.xdim,0,self.Lab.Beam.Lenght))
        plot(Angle,(self.Lab.xdim,0,self.Lab.Beam.Lenght))
        plot(Displacement,(self.Lab.xdim,0,self.Lab.Beam.Lenght))

        
    def __str__(self) -> str:
        print('###### '+ self.LoadCaseName + ' ######')
        print('Has '+ str(len(self.Forces))+ ' force/s:')
        for i in self.Forces:
            print('  -',i)
        print('Has '+ str(len(self.Moments))+ ' Moment/s:')
        for i in self.Moments:
            print('  -',i)
        print('Has '+ str(len(self.FixedSupports))+ ' Fixed Support/s:')
        for i in self.FixedSupports:
            print('  -',i)
        print('Has '+ str(len(self.PinnedSupports))+ ' Pinned Support/s:')
        for i in self.PinnedSupports:
            print('  -',i)
        print('Has '+ str(len(self.VerticalSliderSupports))+ ' Vertical Slider Support/s:')
        for i in self.VerticalSliderSupports:
            print('  -',i)
        return'#################'

class Force():
    def __init__(self,ForceName,Magnitude,Position,LoadCaseEnvironment) -> None:
        self.Name = ForceName
        self.Magnitude = Magnitude
        self.Position = Position
        self.LoadCase = LoadCaseEnvironment
        self.Equation = BuildMacauly(self)

    def __str__(self) -> str:
        return self.Name +': ' + str(self.Magnitude) + '[Nm] at '+ str(self.Position) + '[m]'

class Moment():
    def __init__(self,MomentName,Magnitude,Position,LoadCaseEnvironment) -> None:
        self.Name = MomentName
        self.Magnitude = Magnitude
        self.Position = Position
        self.LoadCase = LoadCaseEnvironment
        self.Equation = BuildMacauly(self)

    def __str__(self) -> str:
        return self.Name +': ' + str(self.Magnitude) + '[N] at '+ str(self.Position) + '[m]'

class ReactionForce():
    def __init__(self,ForceName,Position,LoadCase) -> None:
        self.Name = ForceName
        self.Magnitude = sp.symbols('F_'+self.Name)
        self.Position = Position
        self.LoadCase = LoadCase
        self.Equation = BuildMacauly(self)
        self.LoadCase.UnknownVariables.append(self.Magnitude)

    def __str__(self) -> str:
        return self.Name +': ' + str(self.Magnitude) + '[Nm] at '+ str(self.Position) + '[m]'

class ReactionMoment():
    def __init__(self,MomentName,Position,LoadCase) -> None:
        self.Name = MomentName
        self.Magnitude = sp.symbols('M_'+self.Name)
        self.Position = Position
        self.LoadCase = LoadCase
        self.Equation = BuildMacauly(self)
        self.LoadCase.UnknownVariables.append(self.Magnitude)


    def __str__(self) -> str:
        return self.Name +': ' + str(self.Magnitude) + '[N] at '+ str(self.Position) + '[m]'

class FixedSupport():
    def __init__(self,FixedSupportName,Position,Loadcase) -> None:
        self.Name = FixedSupportName
        self.Loadcase = Loadcase
        self.Force = ReactionForce(FixedSupportName ,Position,self.Loadcase)
        self.Moment = ReactionMoment(FixedSupportName ,Position,self.Loadcase)
        self.Position = Position

        self.Loadcase.Forces.append(self.Force)
        self.Loadcase.Moments.append(self.Moment)
        self.Loadcase.BoundaryConditions.append(DisplacementBoundaryCondition(Position))
        self.Loadcase.BoundaryConditions.append(AngleBoundaryCondition(Position))

    def __str__(self) -> str:
        print(self.Name+':')
        print('    -',self.Force)
        print('    -',self.Moment)
        return '_______________'

class PinnedSupport():
    def __init__(self,PinnedSupportName,Position,Loadcase) -> None:
        self.Name = PinnedSupportName
        self.Loadcase = Loadcase
        self.Force = ReactionForce(PinnedSupportName ,Position,self.Loadcase)
        self.Position = Position

        self.Loadcase.Forces.append(self.Force)
        self.Loadcase.BoundaryConditions.append(DisplacementBoundaryCondition(Position))


    def __str__(self) -> str:
        print(self.Name+':')
        print('    -',self.Force)
        return '_______________'

class VerticalSliderSupport():
    def __init__(self,VerticalSliderSupportName,Position,Loadcase) -> None:
        self.Name = VerticalSliderSupportName
        self.Loadcase = Loadcase
        self.Moment = ReactionMoment(VerticalSliderSupportName,Position,self.Loadcase)
        self.Position = Position

        self.Loadcase.Moments.append(self.Moment)
        self.Loadcase.BoundaryConditions.append(AngleBoundaryCondition(Position))


    def __str__(self) -> str:
        print(self.Name+':')
        print('    -',self.Moment)
        return '_______________'

class ShearForceBoundaryCondition():
    def __init__(self,Position, ShearForce = 0) -> None:
        self.Position = Position
        self.Value = ShearForce

class BendingMomentBoundaryCondition():
    def __init__(self,Position, BendingMoment = 0) -> None:
        self.Position = Position
        self.Value = BendingMoment
        
class AngleBoundaryCondition():
    def __init__(self,Position, AngleDeg = 0) -> None:
        self.Position = Position
        self.Value = radians(AngleDeg)

class DisplacementBoundaryCondition():
    def __init__(self,Position, Displacement = 0) -> None:
        self.Position = Position
        self.Value = Displacement

def BuildMacauly(ForceObject):
    BeamLenght = ForceObject.LoadCase.Lab.Beam.Lenght
    eq = 0
    if ForceObject.Position == BeamLenght:
        eq = sp.Piecewise((0, ForceObject.LoadCase.Lab.xdim<ForceObject.Position),(ForceObject.Magnitude, ForceObject.LoadCase.Lab.xdim>=ForceObject.Position))
    else:
        eq = sp.Piecewise((0, ForceObject.LoadCase.Lab.xdim<=ForceObject.Position),(ForceObject.Magnitude, ForceObject.LoadCase.Lab.xdim>ForceObject.Position))
    return eq
