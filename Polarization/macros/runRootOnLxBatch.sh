#!/bin/sh

export STAGE_HOST=castorcms
export STAGE_SVCCLASS=cmst3

castorDir=/castor/cern.ch/user/h/hwoehri/cmst3/Polarization/Acceptance/JpsiGun_WithFSR/5Dec2010/Histos
cmsswDir=/afs/cern.ch/user/h/hwoehri/scratch0/Polarization/Acceptances/CMSSW_3_8_6/src
macroDir=/afs/cern.ch/user/h/hwoehri/scratch0/Polarization/hWoehri/Polarization/
macroName=runGeomAcc.C
 
for ((job=0;job<1;job++));
  do
  echo "JOB "${job}
  name=processTTree_${job}
#Start to write the script
cat > job_${name}.sh << EOF
#!/bin/sh
cd $cmsswDir
eval \`scramv1 runtime -sh\`
cd -

workDir=$PWD
cp ${macroDir}/macros/runGeomAcc.C .
cp ${macroDir}/macros/GeomAcc.C .
cp ${macroDir}/macros/GeomAcc.h .
cp ${macroDir}/macros/calcPol.C .
cp ${macroDir}/interface/rootIncludes.inc .
cp ${macroDir}/interface/commonVar.h .

echo "before running"
ls -l
root -b -q '${macroName}+(${job})'
echo "after running"
ls -l

rfcp *.root ${castorDir}/

rm -f *.root
EOF
chmod 755 job_${name}.sh
bsub -q 1nh -J $name $PWD/job_${name}.sh

done

