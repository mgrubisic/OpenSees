# home/neeshub/mgrubisic/data/sessions/76987/portal_frame3.tcl

wipe
wipeReliability

# +---------------------+
# |     MODEL BUILD     |
# +---------------------+

model basic -ndm 2 -ndf 3

set dataDir Output;	# name of output folder
file mkdir $dataDir;

node 1 0    0
node 2 0    144
node 3 240  144
node 4 240  0

fix 1   1 1 1
fix 4   1 1 1

set E   30000
set Ag  25
set Ig  1500
set Ac  29
set Ic  2000

section Elastic 1 $E $Ag $Ig
section Elastic 2 $E $Ac $Ic

geomTransf Linear 1

element forceBeamColumn 1 1 2 1 Lobatto 2 3
element forceBeamColumn 2 2 3 1 Lobatto 1 3
element forceBeamColumn 3 3 4 1 Lobatto 2 3

set P 25.0
set w 1.0e-1
pattern Plain 1 Linear {
    load 2 $P 0 0
    eleLoad -ele 2 -type beamUniform -$w 
}

# recorder Node -file $dataDir/portal3.out -time -node 2 -dof 1 disp

# ======== Static Analysis Setup ========
system          BandGeneral
constraints     Transformation
numberer        RCM

set Tol         1.e-8;              
set maxNumIter  6;  
set printFlag   0;            
set TestType    EnergyIncr;	
test $TestType $Tol $maxNumIter $printFlag;

algorithm       Newton
integrator      LoadControl 1
# integrator      LoadControl 0.1
analysis        Static 
# analyze 10

# return

# +------------------------------+
# |     RELIABILITY ANALYSIS     |
# +------------------------------+

reliability

randomVariable 11 lognormal -mean $E -stdv [expr 0.1*$E]
randomVariable 22 normal -mean $P -stdv [expr 0.2*$P]
randomVariable 33 normal -mean 0 -stdv 1
randomVariable 44 normal -mean -$w -stdv [expr abs(0.2*$w)]

parameter 11 randomVariable 11 element 1 E
addToParameter 11 element 2 E
addToParameter 11 element 3 E
parameter 22 randomVariable 22 loadPattern 1 loadAtNode 2 1
parameter 33 randomVariable 33 node 1 coord 1
parameter 44 randomVariable 44 loadPattern 1 elementLoad 2 wy

parameter 55 node 2 disp 1

performanceFunction 76 "0.15-\$par(55)"

# ======== Reliability Analysis Setup ========
sensitivityIntegrator   -static
sensitivityAlgorithm    -computeAtEachStep

randomNumberGenerator           CStdLib
probabilityTransformation       Nataf       -print  3
reliabilityConvergenceCheck     Standard    -e1     1.0e-2  -e2 1.0e-2  -print 1
functionEvaluator               Tcl         -file "analyze 1"
gradientEvaluator               Implicit
searchDirection                 iHLRF
meritFunctionCheck              AdkZhang    -multi  2.0     -add 50     -factor 0.5  
# stepSizeRule                    Armijo      -maxNum 10      -base 0.5   -initial 0.3 5
stepSizeRule                    Fixed       -stepSize 1.0
startPoint                      Mean
findDesignPoint                 StepSearch  -maxNumIter 30; # -printDesignPointX designPointX.out

# ======== Run FORM Analysis ========
runFORMAnalysis  FORMportalframe.out

foreach perf [getLSFTags] {
    puts "Performance Function $perf"
    puts "FORM beta = [format %.7f $betaFORM($perf)]"
    foreach rv [getRVTags] {
   puts "\t x*($rv) = [format %7.4f $designPointXFORM($perf,$rv)], alpha($rv) = [format %.7f $alphaFORM($perf,$rv)], gamma($rv) = [format %.7f $gammaFORM($perf,$rv)]"
    }
}

puts ""
set beta $betaFORM(76);  # betaFORM (lsfTag) gives the beta value for LSF $lsfTag
puts "FORM beta = [format %.7f $betaFORM($perf)]"
puts "FORM Pf = [getStdNormalCDF -$beta]";     # pf = Phi(-beta)


# Run the SORM analysis
findCurvatures		firstPrincipal
runSORMAnalysis		SORMportalframe1.out

hessianEvaluator			FiniteDifference -pert 1000
probabilityTransformation   Nataf           -print 0
findCurvatures				curvatureFitting
runSORMAnalysis				SORMportalframe2.out


return


# +--------------------------------+
# |     MONTE CARLO SIMULATION     |
# +--------------------------------+

set Ntrials 100
set Nfail 0

set startTime [clock clicks -milliseconds]

for {set i 1} {$i <= $Ntrials} {incr i} {  
    reset
    set U ""
    foreach rv [getRVTags] {
        set val [expr rand()]
        lappend U [getStdNormalInverseCDF $val]
        }
        
        set X [transformUtoX $U]
        set irv 0
        
        foreach rv [getRVTags] {
        updateParameter $rv [lindex $X $irv]
        incr irv
        }
        
        recorder Node -file $dataDir/portal3$i.out -time -node 2 -dof 1 disp
        
        integrator      LoadControl 0.1
        analyze 10

        # print -node -flag 1 2
        set D "\[nodeDisp 2 1\]"
        if {[expr 0.15 - $D] <= 0.0} {incr Nfail}
    }
    

set MCSPf [expr double($Nfail)/($Ntrials)];
puts ""
puts "Number of Trials = $Ntrials"
puts "Monte Carlo Pf = [format %.10f  $MCSPf]"
puts ""
puts "FORM beta = [format %.7f $betaFORM($perf)]"
puts "FORM Pf = [getStdNormalCDF -$beta]";     # pf = Phi(-beta)

puts ""
puts "========== End of Reliability Analysis =========="

set finishTime   [clock clicks -milliseconds];
set timeSeconds  [expr ($finishTime-$startTime)/1000];
set timeMinutes  [expr ($timeSeconds/60)];
set timeHours    [expr ($timeSeconds/3600)];
set timeMinutes  [expr ($timeMinutes - $timeHours*60)];
set timeSeconds  [expr ($timeSeconds - $timeMinutes*60 - $timeHours*3600)];

puts "\n--------------------------------------------------------";
#puts "\a";
puts "TOTAL TIME TAKEN: $timeHours hours: $timeMinutes minutes: $timeSeconds seconds";

puts ""
