#!/bin/csh

########################
##  calculation setup ##
########################
#echo "set memory quantity"
#set mem = $<
#echo ""

echo "set the DFT functional"
set dft = $<
echo ""

#echo "set gen and pseudo"
#set ecp = $<
#echo ""

#echo "set charge"
#set charge = $<
#echo ""

#######################
## standard keywords ##
#######################

## SP ##
#set keywords="#p gfinput $dft/gen pseudo=read nosymm int=ultrafine"
set keywords="#p gfinput wb97xd/def2svp nosymm int=ultrafine" #quick sp


## SP + Freq ##
#set keywords="#p freq=(noraman,HPModes) $dft/gen pseudo=read nosymm int=ultrafine"
#set keywords="#p empiricaldispersion=gd3bj freq=(noraman,HPModes) $dft/gen pseudo=read nosymm int=ultrafine"

## Opt ##
#set keywords="#p opt empiricaldispersion=gd3bj $dft/gen pseudo=read nosymm int=ultrafine"
#set keywords="#p opt $dft/gen pseudo=read nosymm int=ultrafine"

## Opt + Freq ##
#set keywords="#p opt freq=(noraman,HPModes) $dft/gen pseudo=read nosymm int=ultrafine"
#set keywords="#p opt empiricaldispersion=gd3 freq=(noraman,HPModes) $dft/def2svp nosymm int=ultrafine" #quick opt
#set keywords="#p opt freq=(noraman,HPModes) $dft/def2svp nosymm int=ultrafine" #quick opt

## OptTS + Freq ##
#set keywords="#p opt=(ts,calcfc,noeigentest, maxstep=5) empiricaldispersion=gd3bj freq=(noraman,HPModes) $dft/gen pseudo=read nosymm int=ultrafine"
#set keywords="#p opt=(ts,calcfc,noeigentest, maxstep=5) freq=(noraman,HPModes) $dft/gen pseudo=read nosymm int=ultrafine"
#set keywords="#p opt=(ts,calcfc,noeigentest, maxstep=5) empiricaldispersion=gd3 freq=(noraman,HPModes) $dft/def2svp nosymm int=ultrafine" #quickts
#set keywords="#p opt=(ts,calcfc,noeigentest, maxstep=5) freq=(noraman,HPModes) $dft/def2svp nosymm int=ultrafine" #quickts

## OptRedundant ##
#set keywords="#p opt=modredundant empiricaldispersion=gd3bj $dft/gen pseudo=read nosymm int=ultrafine"

########################
########################

foreach JOB ( $argv ) 
	## link 0 ##
	echo "%nprocs=8" > $JOB:r.com
	echo "%mem=15gb" >> $JOB:r.com
#	echo "%chk=$JOB:r" >> $JOB:r.com
	## keywords ##
	echo $keywords >> $JOB:r.com
	echo '' >> $JOB:r.com
	## title ##
	echo $JOB:r >> $JOB:r.com
	echo '' >> $JOB:r.com
	## molecule specification ##
	echo "0 1" >> $JOB:r.com
	cat $JOB | sed '1,2 d' >> $JOB:r.com

	## gen and pseudo specification ##
#	echo '' >> $JOB:r.com
#	echo "C H O N 0" >> $JOB:r.com
#	echo 'def2svp' >> $JOB:r.com
#	echo '****' >> $JOB:r.com
#	echo 'Mo 0' >>  $JOB:r.com
#	echo $ecp >>  $JOB:r.com
#	echo '****' >>  $JOB:r.com
#	echo '' >>  $JOB:r.com
#	echo 'Mo 0' >>  $JOB:r.com
#	echo $ecp >> $JOB:r.com
	
	echo '' >>  $JOB:r.com
	echo '' >>  $JOB:r.com
end

echo "The inputs have been created"
echo "Check if every center has a basis-set defined"
 
exit
