#!/usr/bin/wish

###################################################################################################
#    Copyright 2006-2017, S. Merkel, Universite Lille 1, France
#    Contact: email sebastien.merkel at univ-lille1.fr
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
###################################################################################################

###################################################################################################
#
# Requirements:
#	- tcl/tk: this code is in tcl/tk
#	- gnuplot: called for plotting
#	- convert: called to convert EPS to GIF images
#
# Assumptions:
#	assumes that the *.dif input file has the same orientations for all peaks
#		OK: 100 at 90 and 0 degrees, 002 at 90 and 0 degrees, 101 at 90 and 0 degrees
#		NOT OK: 100 at 90 and 0 degrees, 002 at 90 and 0 degrees, 101 at 90 and 45 degrees
#
# Global variables
#	- nsteps: number of steps
#	- stress: array with stresses
#		stress($i,11): sigma_11 for step i
#		stress($i,22): sigma_22 for step i
#		stress($i,33): sigma_33 for step i
#		stress($i,P): pressure for step i (av. s11, s22, s33)
#		stress($i,t): differential stress for step i (s33-s11)
#	- strain
#		strain($i,11): eps_11 for step i
#		strain($i,22): eps_22 for step i
#		strain($i,33): eps_33 for step i
#		strain($i,epratio): ratio of total strain accomodated by elasticity
#	- nsystems: number of deformation mechanisms
#	- activities: array with activities of deformation mechanisms
#		activities($i,$j): activity of system number $j at step $i
#	- nhkl: number of unique hkl planes studied
#	- hkl: array with hkl indices
#		hkl($i) $hkl for plane $i
#	- nangles: number of unique angle studied
#	- angles: array with psi angles values
#		angles($i) $psi for angle $i
#	- norientations: total number of hkl and orientations measured (same as in *.dif file)
#		should be equal to nhkl*nangles
#	- orientations: array with diffraction information, if the dif file is set properly, 
#     it should be redundant with angles and hkl arrays
#		orientations($i,hkl) hkl indices for orientation $i
#		orientations($i,psi) angle psi for orientation $i
#	- diff: array with simulated lattice strains
#		diff($i,$orientations($j,hkl),$orientations($j,psi)) lattice strain for plane and orientation
#		  defined in the 'orientations' array
#
# Changes
#  06/2008, function to  plot stress (s11, s22, s33) in polycrystal vs eps33
# Changed 04/2008 to use actual linear regressions to calculate Q and eps hydrostatic
# Changed 2010 for EPSC4
# Changed 03/2011: calculation of errors on lattice strains
# Changed 07/2014: adapted for EVPSC
# Changed 23/09/2015 S. Merkel to account for elastic contribution to total activity in plot for system activities
# Changed 02/2017: fixed small bug while parsing psi angle values for a non-hexagonal system
###################################################################################################

##################################################
# linear'regression
#	Object: performs linear regression
#	intput: list of values, each item is a two item list {x y}
#	output: returns a two item list {a b} where y = a x + b
# History
#	Created 4/8/2008 S. Merkel
#	Changed 03/2010 to calculate errors, followed numerical recipies
#	Not so useful since we do not really have errors on calculated d-spacings
# Comment:
# Found at http://wiki.tcl.tk/18055
##################################################
proc linear'regression xys {
    set xsum 0.0; set ysum 0.0
    set xsum2 0.0; set xysum 0.0
    set s 0.
    foreach xy $xys {
        foreach {x y} $xy break
		set s [expr $s+1.]
        set xsum [expr {$xsum + $x}]
        set ysum [expr {$ysum + $y}]
        set xsum2 [expr {$xsum2 + $x*$x}]
        set xysum [expr {$xysum + $x*$y}]
    }
    set delta [expr $s*$xsum2-$xsum*$xsum]
    set a [expr ($xsum2*$ysum-$xsum*$xysum)/$delta]
    set b [expr ($s*$xysum-$xsum*$ysum)/$delta]
	# If only 2 data points, errors do not make any sense
	if {[llength $xys] > 2} {
		set chi2 0.
		foreach xy $xys {
			foreach {x y} $xy break
			set t [expr $y-$a-$x*$b]
			set chi2 [expr {$chi2 + $t*$t}]
		}
		set sa [expr sqrt($xsum2/$delta)*sqrt($chi2/($s-2.))]
		set sb [expr sqrt($s/$delta)*sqrt($chi2/($s-2.))]
	} else {
		set sa 0
		set sb 0
	}
    list $a $b $sa $sb
    
    
#     set xm [expr {$xsum/[llength $xys]}]
#     set ym [expr {$ysum/[llength $xys]}]
#     set xsum 0.0; set ysum 0.0
#     foreach xy $xys {
#         foreach {x y} $xy break
#         set dx [expr {$x - $xm}]
#         set dy [expr {$y - $ym}]
#         set xsum [expr {$xsum + $dx * $dy}]
#         set ysum [expr {$ysum + $dx * $dx}]
#     }
#     set b [expr {$xsum / $ysum}]
#     set a [expr {$ym - $b * $xm}]
#     list $a $b
 }

##################################################
# dpQ
#	Object: fits the relation eps = eps0 + (1+eps0) (Q (1 - 3 cos^2 psi))
#	intput: list of values, each item is a two item list {psi eps}
#	output: returns a two item list {eps0 Q}
# History
#	Created 4/8/2008 S. Merkel
#	Changed 03/2011 S. Merkel to return errors
# Comment:
##################################################
proc dpQ psids {
	# puts $psids
	set pi [expr acos(-1.)]
	set xys {}
	foreach psid $psids {
		foreach {psi d} $psid break;
		set cospsi [expr cos($psi*$pi/180.)]
		set x [expr 1.-3.*$cospsi*$cospsi]
		lappend xys [list $x $d]
    }
	# puts $xys
    foreach {a b sa sb} [linear'regression $xys] break
	set dp $a
	set Q [expr $b/($a+1.)]
	set t1 [expr $sb/($a+1.)]
	set t2 [expr $b*$sa/(($a+1.)*($a+1.))]
	set dQ [expr sqrt($t1*$t1+$t2*$t2)]
	list $dp $Q $sa $dQ
}

###################################################################################################
#
# Input procedudes: procedures used to parse the EPSC output
#
###################################################################################################

##################################################
# lineToArray
#	Object: converts a line of space-separated data to an array
#	intput: line of data
#	output: array of data
# History
#	Created 10/10/2006 S. Merkel
##################################################
proc lineToArray {line} {
	regsub -all "D" $line "e" newline
	regsub -all "\"" $newline "" newline
	regsub -all "\'" $newline "" newline
	set dataarray [concat $newline]
	return $dataarray
}

##################################################
# countsteps
#	Object: counts the number of steps in the calculation (using the activity file)
#	intput:
#	output: sets the nsteps global variable
# History
#	Created 10/10/2006 S. Merkel
#   Changed 04/12/2008 S. Merkel to adapt to EPSC 4 (columns are different)
#   Changed 07/11/2014 S. Merkel to adapt to EVPSC
#   Changed 23/09/2015 S. Merkel, STR_STR has one more component
##################################################
proc countsteps {} {
	global nsteps
	set fichier [open "STR_STR.OUT" "r"]
	# label lines at the top
	gets $fichier data
	# count the lines
	set nsteps  0
	while {[eof $fichier] == 0} {
		gets $fichier data
		set dataarray [lineToArray $data]
		set ncols [llength $dataarray]
		if {$ncols == 15} {
			incr nsteps
		}
	}
	close $fichier
	puts "Parsed STR_STR.OUT. Found number of steps: $nsteps."
}

##################################################
# readstresses
#	Object: parses the stress and strain file (step, strains, stresses)
#	intput:
#	output: sets the stress and strain global variables
# History
#	Created 10/10/2006 S. Merkel
#   Changed 04/12/2008 S. Merkel to adapt to EPSC 4 (columns are different)
#   Changed 07/11/2014 S. Merkel to adapt to EVPSC
#   Changed 23/09/2015 S. Merkel, reads elastic vs. plastic deformation component
##################################################
proc readstresses {} {
	global nsteps
	global stress
	global strain
	set fichier [open "STR_STR.OUT" "r"]
	# 6 label lines at the top
	gets $fichier data
	# load data
	for {set i 0} {$i < $nsteps} {incr i} {
		gets $fichier data
		set dataarray [lineToArray $data]
		set strain($i,11) [expr -[lindex $dataarray 2]]
		set strain($i,22) [expr -[lindex $dataarray 3]]
		set strain($i,33) [expr -[lindex $dataarray 4]]
		set stress($i,11) [expr -[lindex $dataarray 8]]
		set stress($i,22) [expr -[lindex $dataarray 9]]
		set stress($i,33) [expr -[lindex $dataarray 10]]
		set stress($i,P) [expr ($stress($i,11)+$stress($i,22)+$stress($i,33))/3.]
		set stress($i,t) [expr $stress($i,33)-$stress($i,11)]
		set strain($i,epratio) [expr [lindex $dataarray 14]]
	}
	close $fichier
	puts "Parsed STR_STR.OUT for stresses."
}

##################################################
# readactivities
#	Object: parses the activity files file (slip system activites)
#	intput: phase: phase number (1, 2, 3...)
#	output: sets the activities and nsystems global variable
# History
#	Created 10/10/2006 S. Merkel
#   Changed 12/04/2008 S. Merkel to read EPSCV4 format
#   Changed 07/11/2014 S. Merkel to adapt to EVPSC
##################################################
proc readactivities {phase} {
	global nsteps
	global activities
	global nsystems
	set fichier [open "ACT_PH${phase}.OUT" "r"]
	gets $fichier data
	set dataarray [lineToArray $data]
	set n [llength $dataarray]
	set nsystems 0
	for {set i 0} {$i < $n} {incr i} {
		set txt [lindex $dataarray $i]
		if {[string match "*MO*" $txt]} {
			if {$nsystems == 0} {set mode1 $i}
			incr nsystems
			
		}
	}
	# activities at step 0 are null
	for {set j 0} {$j < $nsystems} {incr j} {
		set activities(0,$j) 0
	}
	for {set i 1} {$i <= $nsteps-1} {incr i} {
		gets $fichier data
		set dataarray [lineToArray $data]
		for {set j 0} {$j < $nsystems} {incr j} {
			set activities($i,$j) [lindex $dataarray [expr $j+$mode1]]
		}
	}
	close $fichier
	puts "Parsed ACT_PH${phase}.OUT for slip system activities."
}

##################################################
# readdiffraction
#	Object: parses the diffraction file (diffraction data)
#	intput: diffraction configuration file, phase number (1, 2, 3...)
#	output: sets the norientations, nstepOrientations, orientations, and diff global variables
# History
#	Created 10/10/2006 S. Merkel
# Comment:
#	This routine depends on simulation! It uses the .dif file 
#	to set a figure out what are the observations
# Changed
#   Changed 07/11/2014 S. Merkel to adapt to EVPSC
#  04/12/2008 to accomodate for V4 EPSC
#   Changed 07/11/2014 S. Merkel to adapt to EVPSC
##################################################
proc readdiffraction {diffile phase} {
	global norientations
	global orientations
	global nhkl
	global hkl
	global angles
	global nangles
	global diff
	global nsteps
	# Reading cobal.dif to extract orientation informations
	set fichier [open $diffile "r"]
	gets $fichier data
	gets $fichier data
	gets $fichier data
	set dataarray [lineToArray $data]
	set norientations [lindex $dataarray 0]
	set hcpornot [lindex $dataarray 1]
	if {$hcpornot == 2} {
		set lcol 3
		set psicol 4
	} else {
		set lcol 2
		set psicol 3
	}
	gets $fichier data
	gets $fichier data
	set nhkl 0
	set nangles 0
	for {set i 0} {$i < $norientations} {incr i} {
		gets $fichier data
		set dataarray [lineToArray $data]
		set h [lindex $dataarray 0]
		set k [lindex $dataarray 1]
		set l [lindex $dataarray $lcol]
		set psi [expr int([lindex $dataarray $psicol])]
		set thishkl "$h$k$l"
		set orientations($i,hkl) $thishkl
		set orientations($i,psi) $psi
		set found 0
		for {set j 0} {$j < $nhkl} {incr j} {
			if {$hkl($j) == $thishkl} { set found 1 }
		}
		if {$found == 0} {
			set hkl($nhkl) $thishkl
			set nhkl [expr $nhkl+1]
		}
		set found 0
		for {set j 0} {$j < $nangles} {incr j} {
			if {$angles($j) == $psi} { set found 1 }
		}
		if {$found == 0} {
			set angles($nangles) $psi
			set nangles [expr $nangles+1]
		}
	}
	puts "Read $diffile. Found $nhkl hkl and $nangles angles"
	# for {set j 0} {$j < $nhkl} {incr j} { puts "  $hkl($j)" }
	for {set j 0} {$j < $nangles} {incr j} { puts "  $angles($j)" }
	close $fichier
	# reading diffraction data
	set fichier [open "LAT_STR${phase}.OUT" "r"]
	for {set i 0} {$i < $nsteps} {incr i} {
		#puts "Looking at step $i"
		# Pulling info on volume fraction
		# puts "Reading step $i"
		gets $fichier data
		set dataarray [lineToArray $data]
		for {set j 0} {$j < $norientations} {incr j} {
			set ttt [lindex $dataarray [expr $j+7]]
			# puts "Trying to work with $ttt. Step is $i, diff is $j"
			# There used to be a division by 1e6 for EPSC, but it does not seem necessary here
			set diff($i,$orientations($j,hkl),$orientations($j,psi)) [expr [lindex $dataarray [expr $j+7]]/1.e0]
			 #puts "Reading data for $i $orientations($j,hkl) $orientations($j,psi): $diff($i,$orientations($j,hkl),$orientations($j,psi)) "
		}
	}
	close $fichier
	puts "Parsed LAT_STR${phase}.OUT for lattice strains."
}

###################################################################################################
#
# Gnuplot interface: various subroutines used to create plots with gnuplot (window, eps, gif...)
#
###################################################################################################

##################################################
# gnuplotPostscript
#	Object: uses gnuplot to create a B&W postscript plot
#	intput:
#		- indexplot: index of the gnuplot input (array inputGnuplot)
#		- filename: root for the output file (.eps will be added automatically)
#	output:
# History
#	Created 10/11/2006 S. Merkel
# Comment:
##################################################
proc gnuplotPostscript {indexplot filename} {
	global inputGnuplot
	set gp [open "|gnuplot" r+]
	puts $gp "set term postscript eps enhanced lw 2 \"Helvetica\" 18"
	puts $gp "set output \"$filename.eps\""
	puts $gp $inputGnuplot($indexplot)
	close $gp
}

##################################################
# gnuplotPostscriptColor
#	Object: uses gnuplot to create a color postscript plot
#	intput:
#		- indexplot: index of the gnuplot input (array inputGnuplot)
#		- filename: root for the output file (.eps will be added automatically)
#	output:
# History
#	Created 10/11/2006 S. Merkel
# Comment:
##################################################
proc gnuplotPostscriptColor {indexplot filename} {
	global inputGnuplot
	set gp [open "|gnuplot" r+]
	puts $gp "set term postscript eps enhanced color lw 2 \"Helvetica\" 18"
	puts $gp "set output \"$filename.eps\""
	puts $gp $inputGnuplot($indexplot)
	close $gp
}

##################################################
# gnuplotGif
#	Object: uses gnuplot and convert to create a GIF plot
#	intput:
#		- indexplot: index of the gnuplot input (array inputGnuplot)
#		- filename: root for the output file (.gif will be added automatically)
#	output:
# History
#	Created 10/11/2006 S. Merkel
# Comment:
##################################################
proc gnuplotGif {indexplot filename} {
	global inputGnuplot
	set gp [open "|gnuplot" r+]
	puts $gp "set term postscript eps enhanced color lw 2 \"Helvetica\" 18"
	puts $gp "set output \"tmp.eps\""
	puts $gp $inputGnuplot($indexplot)
	close $gp
	exec convert  -density 144 tmp.eps "$filename.gif"
	file delete "tmp.eps"
}

##################################################
# clearGnuplot 
#	Object: gets rid of a plot window and clear gnuplot input cache
#	intput:
#		- indexplot: index of the plot window
#	output:
# History
#	Created 10/11/2006 S. Merkel
# Comment:
##################################################
proc clearGnuplot {indexPlot} {
	global inputGnuplot
	unset inputGnuplot($indexPlot)
	destroy .c${indexPlot}
}

##################################################
# startGnuplot
#	Object: creates a plot window with plot and export buttons
#	intput:
#		- input: gnuplot commands
#		- windowtitle
#		- filename: used to choose filename for GIF or EPS output
#	output:
# History
#	Created 10/11/2006 S. Merkel
# Comment:
##################################################
proc startGnuplot {input windowtitle filename} {
	global indexPlot
	global inputGnuplot
	# Output to /dev/null is important to avoid garbage (when curve fitting for instance)
	# Have a look at http://perso.numericable.fr/~mipicard/sorties.html to adapt to windows
	set gp [open "|gnuplot >& /dev/null" w]
	puts $gp "set term tk"
	puts $gp "set output 'resultat.tk'"
	puts $gp $input
	close $gp
	incr indexPlot
	set inputGnuplot($indexPlot) $input
	toplevel .c${indexPlot}
    wm title .c${indexPlot} $windowtitle
    #wm minsize .c${indexPlot} 500 400
	canvas .c${indexPlot}.plot -width 500 -height 400 -background "#FFFFFF"
	frame .c${indexPlot}.buttons
	button .c${indexPlot}.buttons.exportPS -text "Postscript (NB)"  -width 15 \
			-command "gnuplotPostscript \"$indexPlot\" \"$filename\""
	button .c${indexPlot}.buttons.exportPSColor -text "Postscript (Color)"  -width 15 \
			-command "gnuplotPostscriptColor \"$indexPlot\" \"$filename\""
	button .c${indexPlot}.buttons.exportGif -text "Gif"  -width 15 \
			-command "gnuplotGif \"$indexPlot\" \"$filename\""
	button .c${indexPlot}.buttons.close -text "Close"  -width 15 \
			-command "clearGnuplot \"${indexPlot}\""
	pack .c${indexPlot}.plot
    pack .c${indexPlot}.buttons.close -padx 5 -pady 5 -side left
    pack .c${indexPlot}.buttons.exportPS -padx 5 -pady 5 -side left
    pack .c${indexPlot}.buttons.exportPSColor -padx 5 -pady 5 -side left
    pack .c${indexPlot}.buttons.exportGif -padx 5 -pady 5 -side left
    pack .c${indexPlot}.buttons
	source resultat.tk
	gnuplot .c${indexPlot}.plot
}

###################################################################################################
#
# Plotting procedudes: procedures used to plot results of the EPSC calculations
#
###################################################################################################


##################################################
# plotStressesvsEps33
#	Object: creates a plot of s11, s22, s33 vs. eps33
#	intput:
#	output:
# History
#	Created 06/06/2008 S. Merkel
# Comment:
##################################################
proc plotStressesvsEps33 {} {
	global nsteps
	global stress
	global strain
	set output [open "stressVsEps33.tmp" "w"]
	for {set i 0} {$i < $nsteps} {incr i} {
		puts $output "$strain($i,33)\t$stress($i,11)\t$stress($i,22)\t$stress($i,33)"
	}
	close $output
	set gnuplot "set ylabel \"Sigma (GPa)\"\nset xlabel \"Eps33\"\nplot \"stressVsEps33.tmp\" using 1:2 title \"s11 vs Eps33\", \"stressVsEps33.tmp\" using 1:3 title \"s22 vs Eps33\", \"stressVsEps33.tmp\" using 1:4 title \"s33 vs Eps33\""
	startGnuplot $gnuplot "Sigma vs Eps33" "stressVsEps33"
}


##################################################
# plotPvsEps33
#	Object: creates a plot of pressure vs. eps33
#	intput:
#	output:
# History
#	Created 10/11/2006 S. Merkel
# Comment:
##################################################
proc plotPvsEps33 {} {
	global nsteps
	global stress
	global strain
	set output [open "pVsEps33.tmp" "w"]
	for {set i 0} {$i < $nsteps} {incr i} {
		puts $output "$strain($i,33)\t$stress($i,P)"
	}
	close $output
	set gnuplot "set ylabel \"P (GPa)\"\nset xlabel \"Eps33\"\nplot \"pVsEps33.tmp\" using 1:2 title \"P vs Eps33\""
	startGnuplot $gnuplot "P vs Eps33" "pVsEps33"
}

##################################################
# plotTvsEps33
#	Object: creates a plot of differential stress vs. eps33
#	intput:
#	output:
# History
#	Created 10/11/2006 S. Merkel
# Comment:
##################################################
proc plotTvsEps33 {} {
	global nsteps
	global stress
	global strain
	set output [open "tVsEps33.tmp" "w"]
	for {set i 0} {$i < $nsteps} {incr i} {
		puts $output "$strain($i,33)\t$stress($i,t)"
	}
	close $output
	set gnuplot "set ylabel \"t (GPa)\"\nset xlabel \"Eps33\"\nplot \"tVsEps33.tmp\" using 1:2 title \"t vs Eps33\""
	startGnuplot $gnuplot "t vs Eps33" "tVsEps33"
}

##################################################
# plotTvsP
#	Object: creates a plot of differential stress vs. pressure
#	intput:
#	output:
# History
#	Created 10/11/2006 S. Merkel
# Comment:
##################################################
proc plotTvsP {} {
	global nsteps
	global stress
	global strain
	set output [open "tVsP.tmp" "w"]
	for {set i 0} {$i < $nsteps} {incr i} {
		puts $output "$stress($i,P)\t$stress($i,t)"
	}
	close $output
	set gnuplot "set ylabel \"t (GPa)\"\nset xlabel \"P (GPa)\"\nplot \"tVsP.tmp\" using 1:2 title \"t vs Eps33\""
	startGnuplot $gnuplot "t vs P" "tVsP"
}

##################################################
# plotAvsEps33
#	Object: creates a plot of slip system activities vs. eps33
#	intput:
#	output:
# History
#	Created 10/11/2006 S. Merkel
#	Changed 23/09/2015 S. Merkel to account for elastic contribution
# Comment:
##################################################
proc plotAvsEps33 {} {
	global nsteps
	global strain
	global activities
	global nsystems
	global strain

	set output [open "actVsEps33.tmp" "w"]
	for {set i 0} {$i < $nsteps} {incr i} {
		set line "$strain($i,33)"
		set line "$line\t$strain($i,epratio)"
		for {set j 0} {$j < $nsystems} {incr j} {
			set act [expr $activities($i,$j)*(1.-$strain($i,epratio))]
			set line "$line\t$act"
		}
		puts $output $line
	}
	close $output
	set gnuplot "set ylabel \"Activity\"\nset xlabel \"Eps33\""
	set gnuplot "$gnuplot\nplot \"actVsEps33.tmp\" using 1:2 title \"Elastic\""
	for {set j 0} {$j < $nsystems} {incr j} {
		set col [expr $j+3]
		set nsys [expr $j+1]
		set gnuplot "$gnuplot, \"actVsEps33.tmp\" using 1:$col title \"System $nsys\""
	}
	puts $gnuplot
	startGnuplot $gnuplot "Slip systems activities vs Eps33" "actVsEps33"
}

##################################################
# plotAvsStep
#	Object: creates a plot of slip system activities vs. step number
#	intput:
#	output:
# History
#	Created 10/11/2006 S. Merkel
#	Changed 23/09/2015 S. Merkel to account for elastic contribution
# Comment:
##################################################
proc plotAvsStep {} {
	global nsteps
	global activities
	global nsystems
	global strain

	set output [open "actVsStep.tmp" "w"]
	for {set i 0} {$i < $nsteps} {incr i} {
		set line "$i"
		set line "$line\t$strain($i,epratio)"
		for {set j 0} {$j < $nsystems} {incr j} {
			set act [expr $activities($i,$j)*(1.-$strain($i,epratio))]
			set line "$line\t$act"
		}
		puts $output $line
	}
	close $output
	set gnuplot "set ylabel \"Activity\"\nset xlabel \"Step number\""
	set gnuplot "$gnuplot\nplot \"actVsStep.tmp\" using 1:2 title \"Elastic\""
	for {set j 0} {$j < $nsystems} {incr j} {
		set col [expr $j+3]
		set nsys [expr $j+1]
		set gnuplot "$gnuplot, \"actVsStep.tmp\" using 1:$col title \"System $nsys\""
	}
	startGnuplot $gnuplot "Slip system activities vs step number" "actVsStep"
}

##################################################
# plotAvsP
#	Object: creates a plot of slip system activities vs. pressure
#	intput:
#	output:
# History
#	Created 10/11/2006 S. Merkel
#	Changed 23/09/2015 S. Merkel to account for elastic contribution
# Comment:
##################################################
proc plotAvsP {} {
	global nsteps
	global stress
	global activities
	global nsystems
	global strain

	set output [open "actVsP.tmp" "w"]
	for {set i 0} {$i < $nsteps} {incr i} {
		set line "$stress($i,P)"
		set line "$line\t$strain($i,epratio)"
		for {set j 0} {$j < $nsystems} {incr j} {
			set act [expr $activities($i,$j)*(1.-$strain($i,epratio))]
			set line "$line\t$act"
		}
		puts $output $line
	}
	close $output
	set gnuplot "set ylabel \"Activity\"\nset xlabel \"P (GPa)\""
	set gnuplot "$gnuplot\nplot \"actVsP.tmp\" using 1:2 title \"Elastic\""
	for {set j 0} {$j < $nsystems} {incr j} {
		set col [expr $j+3]
		set nsys [expr $j+1]
		set gnuplot "$gnuplot, \"actVsP.tmp\" using 1:$col title \"System $nsys\""
	}
	startGnuplot $gnuplot "Slip systems activities vs P" "actVsP"
}


##################################################
# plotErrQvsEps33
#	Object: creates a plot of standard deviations for Q values vs. eps33
#	intput:
#	output:
# History
#	Created Mar 24 2010 S. Merkel
# Comment:
############################################################################################

proc plotErrQvsEps33 {} {
	global strain
	global nsteps
	global orientations
	global diff
	global nhkl
	global hkl
	global nangles
	global angles
	set output [open "SigmaQVsEps33.tmp" "w"]
	set pi [expr acos(-1.)]
	for {set i 0} {$i < $nsteps} {incr i} {
		set line "$strain($i,33)"
		for {set j 0} {$j < $nhkl} {incr j} {
			set psids {}
			for {set k 0} {$k < $nangles} {incr k} {
				set psi $angles($k)
				set d $diff($i,$hkl($j),$angles($k))
				lappend psids [list $psi $d]
			}
			#puts $psids
			foreach {dp Q ddp dQ} [dpQ $psids] break
			set line "$line\t$dQ"
		}
		puts $output $line
	}
	close $output
	set gnuplot "set ylabel \"sigma Q(hkl)\"\nset xlabel \"Eps33\"\nset key left top"
	for {set j 0} {$j < $nhkl} {incr j} {
		set col [expr $j+2]
		set nsys [expr $j+1]
		if {$nsys == 1} {
			set gnuplot "$gnuplot\nplot \"SigmaQVsEps33.tmp\" using 1:$col title \"$hkl($j)\""
		} else {
			set gnuplot "$gnuplot, \"SigmaQVsEps33.tmp\" using 1:$col title \"$hkl($j)\""
		}
	}
	startGnuplot $gnuplot "Sigma(Q) vs Eps33" "SigmaQVsEps33"
}

##################################################
# plotQvsStep
#	Object: creates a plot of Q values vs. step
#	intput:
#	output:
# History
#	Created 10/04/2016 S. Merkel from plotQvsEps33
# Comment:
########

proc plotQvsStep {} {
	global strain
	global nsteps
	global orientations
	global diff
	global nhkl
	global hkl
	global nangles
	global angles
	set output [open "QVsStep.tmp" "w"]
	for {set i 0} {$i < $nsteps} {incr i} {
		set line "$i"
		for {set j 0} {$j < $nhkl} {incr j} {
			set psids {}
			for {set k 0} {$k < $nangles} {incr k} {
				set psi $angles($k)
				set d $diff($i,$hkl($j),$angles($k))
				lappend psids [list $psi $d]
			}
			#puts $psids
			foreach {dp Q} [dpQ $psids] break
			#set eps0 $diff($i,$hkl($j),0)
			#set eps90 $diff($i,$hkl($j),90)
			#set Q [expr -($eps90-$eps0)/($eps90+2*$eps0-3)]
			
			set line "$line\t$Q"
		}
		puts $output $line
	}
	close $output
	set gnuplot "set ylabel \"Q(hkl)\"\nset xlabel \"Step\"\nset key left top"
	for {set j 0} {$j < $nhkl} {incr j} {
		set col [expr $j+2]
		set nsys [expr $j+1]
		if {$nsys == 1} {
			set gnuplot "$gnuplot\nplot \"QVsStep.tmp\" using 1:$col title \"$hkl($j)\""
		} else {
			set gnuplot "$gnuplot, \"QVsStep.tmp\" using 1:$col title \"$hkl($j)\""
		}
	}
	startGnuplot $gnuplot "Lattice strain Q vs Step" "QVsStep"
}


##################################################
# plotQvsEps33
#	Object: creates a plot of Q values vs. eps33
#	intput:
#	output:
# History
#	Created 10/12/2006 S. Merkel
#	Modified 4/8/2008 S. Merkel: real linear fit to get Q values
# Comment:
########

proc plotQvsEps33 {} {
	global strain
	global nsteps
	global orientations
	global diff
	global nhkl
	global hkl
	global nangles
	global angles
	set output [open "QVsEps33.tmp" "w"]
	for {set i 0} {$i < $nsteps} {incr i} {
		set line "$strain($i,33)"
		for {set j 0} {$j < $nhkl} {incr j} {
			set psids {}
			for {set k 0} {$k < $nangles} {incr k} {
				set psi $angles($k)
				set d $diff($i,$hkl($j),$angles($k))
				lappend psids [list $psi $d]
			}
			#puts $psids
			foreach {dp Q} [dpQ $psids] break
			#set eps0 $diff($i,$hkl($j),0)
			#set eps90 $diff($i,$hkl($j),90)
			#set Q [expr -($eps90-$eps0)/($eps90+2*$eps0-3)]
			
			set line "$line\t$Q"
		}
		puts $output $line
	}
	close $output
	set gnuplot "set ylabel \"Q(hkl)\"\nset xlabel \"Eps33\"\nset key left top"
	for {set j 0} {$j < $nhkl} {incr j} {
		set col [expr $j+2]
		set nsys [expr $j+1]
		if {$nsys == 1} {
			set gnuplot "$gnuplot\nplot \"QVsEps33.tmp\" using 1:$col title \"$hkl($j)\""
		} else {
			set gnuplot "$gnuplot, \"QVsEps33.tmp\" using 1:$col title \"$hkl($j)\""
		}
	}
	startGnuplot $gnuplot "Lattice strain Q vs Eps33" "QVsEps33"
}



##################################################
# plotErrQvsP
#	Object: creates a plot of sigma(Q) values vs. pressure
#	intput:
#	output:
# History
#	Created 03/2011 S. Merkel
# Comment:
##################################################
proc plotErrQvsP {} {
	global stress
	global nsteps
	global orientations
	global diff
	global nhkl
	global hkl
	global nangles
	global angles
	set output [open "SigmaQVsP.tmp" "w"]
	for {set i 0} {$i < $nsteps} {incr i} {
		set line "$stress($i,P)"
		for {set j 0} {$j < $nhkl} {incr j} {
			set psids {}
			for {set k 0} {$k < $nangles} {incr k} {
				set psi $angles($k)
				set d $diff($i,$hkl($j),$angles($k))
				lappend psids [list $psi $d]
			}
			foreach {dp Q ddp dQ} [dpQ $psids] break
			set line "$line\t$dQ"
		}
		puts $output $line
	}
	close $output
	set gnuplot "set ylabel \"Sigma(Q)(hkl)\"\nset xlabel \"P (GPa)\"\nset key left top"
	for {set j 0} {$j < $nhkl} {incr j} {
		set col [expr $j+2]
		set nsys [expr $j+1]
		if {$nsys == 1} {
			set gnuplot "$gnuplot\nplot \"SigmaQVsP.tmp\" using 1:$col title \"$hkl($j)\""
		} else {
			set gnuplot "$gnuplot, \"SigmaQVsP.tmp\" using 1:$col title \"$hkl($j)\""
		}
	}
	startGnuplot $gnuplot "Sigma(Q) vs pressure" "QVsP"
}

##################################################
# plotQvsP
#	Object: creates a plot of Q values vs. pressure
#	intput:
#	output:
# History
#	Created 10/12/2006 S. Merkel
#	Modified 4/8/2008 S. Merkel: real linear fit to get Q values
# Comment:
##################################################
proc plotQvsP {} {
	global stress
	global nsteps
	global orientations
	global diff
	global nhkl
	global hkl
	global nangles
	global angles
	set output [open "QVsP.tmp" "w"]
	for {set i 0} {$i < $nsteps} {incr i} {
		set line "$stress($i,P)"
		for {set j 0} {$j < $nhkl} {incr j} {
			set psids {}
			for {set k 0} {$k < $nangles} {incr k} {
				set psi $angles($k)
				set d $diff($i,$hkl($j),$angles($k))
				lappend psids [list $psi $d]
			}
			#puts $psids
			foreach {dp Q} [dpQ $psids] break
			#set eps0 $diff($i,$hkl($j),0)
			#set eps90 $diff($i,$hkl($j),90)
			#set Q [expr -($eps90-$eps0)/($eps90+2*$eps0-3)]
			set line "$line\t$Q"
		}
		puts $output $line
	}
	close $output
	set gnuplot "set ylabel \"Q(hkl)\"\nset xlabel \"P (GPa)\"\nset key left top"
	for {set j 0} {$j < $nhkl} {incr j} {
		set col [expr $j+2]
		set nsys [expr $j+1]
		if {$nsys == 1} {
			set gnuplot "$gnuplot\nplot \"QVsP.tmp\" using 1:$col title \"$hkl($j)\""
		} else {
			set gnuplot "$gnuplot, \"QVsP.tmp\" using 1:$col title \"$hkl($j)\""
		}
	}
	startGnuplot $gnuplot "Lattice strain Q vs pressure" "QVsP"
}

##################################################
# plotQvsT
#	Object: creates a plot of Q values vs. differential stress
#	intput:
#	output:
# History
#	Created 10/12/2006 S. Merkel
#	Modified 4/8/2008 S. Merkel: real linear fit to get Q values
# Comment:
##################################################
proc plotQvsT {} {
	global stress
	global nsteps
	global orientations
	global diff
	global nhkl
	global hkl
	global nangles
	global angles
	set output [open "QVsT.tmp" "w"]
	for {set i 0} {$i < $nsteps} {incr i} {
		set line "$stress($i,t)"
		for {set j 0} {$j < $nhkl} {incr j} {
			set psids {}
			for {set k 0} {$k < $nangles} {incr k} {
				set psi $angles($k)
				set d $diff($i,$hkl($j),$angles($k))
				lappend psids [list $psi $d]
			}
			#puts $psids
			foreach {dp Q} [dpQ $psids] break
			#set eps0 $diff($i,$hkl($j),0)
			#set eps90 $diff($i,$hkl($j),90)
			#set Q [expr -($eps90-$eps0)/($eps90+2*$eps0-3)]
			set line "$line\t$Q"
		}
		puts $output $line
	}
	close $output
	set gnuplot "set ylabel \"Q(hkl)\"\nset xlabel \"t (GPa)\"\nset key left top"
	for {set j 0} {$j < $nhkl} {incr j} {
		set col [expr $j+2]
		set nsys [expr $j+1]
		if {$nsys == 1} {
			set gnuplot "$gnuplot\nplot \"QVsT.tmp\" using 1:$col title \"$hkl($j)\""
		} else {
			set gnuplot "$gnuplot, \"QVsT.tmp\" using 1:$col title \"$hkl($j)\""
		}
	}
	startGnuplot $gnuplot "Lattice strain Q vs differential stress" "QVsT"
}

##################################################
# plotEpsHydrovsEps33
#	Object: creates a plot of "hydrostatic strain" values for each plane vs. eps33
#	intput:
#	output:
# History
#	Created 14/11/2007 S. Merkel
#	Modified 4/8/2008 S. Merkel: real linear fit to get hydrostatic strain
# Comment:
##################################################
proc plotEpsHydrovsEps33 {} {
	global strain
	global nsteps
	global orientations
	global diff
	global nhkl
	global hkl
	global nangles
	global angles
	set output [open "EpsHydroVsEps33.tmp" "w"]
	for {set i 0} {$i < $nsteps} {incr i} {
		set line "$strain($i,33)"
		for {set j 0} {$j < $nhkl} {incr j} {
			set psids {}
			for {set k 0} {$k < $nangles} {incr k} {
				set psi $angles($k)
				set d $diff($i,$hkl($j),$angles($k))
				lappend psids [list $psi $d]
			}
			#puts $psids
			foreach {dp Q} [dpQ $psids] break
			#set eps0 $diff($i,$hkl($j),0)
			#set eps90 $diff($i,$hkl($j),90)
			#set epsHydro [expr (2.*$eps90+$eps0)/3.]
			set line "$line\t$dp"
		}
		puts $output $line
	}
	close $output
	set gnuplot "set ylabel \"Hydrostatic strain(hkl)\"\nset xlabel \"Eps33\"\nset key right top"
	for {set j 0} {$j < $nhkl} {incr j} {
		set col [expr $j+2]
		set nsys [expr $j+1]
		if {$nsys == 1} {
			set gnuplot "$gnuplot\nplot \"EpsHydroVsEps33.tmp\" using 1:$col title \"$hkl($j)\""
		} else {
			set gnuplot "$gnuplot, \"EpsHydroVsEps33.tmp\" using 1:$col title \"$hkl($j)\""
		}
	}
	startGnuplot $gnuplot "Hydrostatic strain vs Eps33" "HydroStrainVsEps33"
}

##################################################
# plotEpsHydrovsP
#	Object: creates a plot of "hydrostatic strain" values for each plane  vs. pressure
#	intput:
#	output:
# History
#	Created 14/11/2007 S. Merkel
#	Modified 4/8/2008 S. Merkel: real linear fit to get hydrostatic strain
# Comment:
##################################################
proc plotEpsHydrovsP {} {
	global stress
	global nsteps
	global orientations
	global diff
	global nhkl
	global hkl
	global nangles
	global angles
	set output [open "EpsHydroVsP.tmp" "w"]
	for {set i 0} {$i < $nsteps} {incr i} {
		set line "$stress($i,P)"
		for {set j 0} {$j < $nhkl} {incr j} {
			set psids {}
			for {set k 0} {$k < $nangles} {incr k} {
				set psi $angles($k)
				set d $diff($i,$hkl($j),$angles($k))
				lappend psids [list $psi $d]
			}
			#puts $psids
			foreach {dp Q} [dpQ $psids] break
			#set eps0 $diff($i,$hkl($j),0)
			#set eps90 $diff($i,$hkl($j),90)
			#set epsHydro [expr (2.*$eps90+$eps0)/3.]
			set line "$line\t$dp"
		}
		puts $output $line
	}
	close $output
	set gnuplot "set ylabel \"Hydrostatic strain (hkl)\"\nset xlabel \"P (GPa)\"\nset key right top"
	for {set j 0} {$j < $nhkl} {incr j} {
		set col [expr $j+2]
		set nsys [expr $j+1]
		if {$nsys == 1} {
			set gnuplot "$gnuplot\nplot \"EpsHydroVsP.tmp\" using 1:$col title \"$hkl($j)\""
		} else {
			set gnuplot "$gnuplot, \"EpsHydroVsP.tmp\" using 1:$col title \"$hkl($j)\""
		}
	}
	startGnuplot $gnuplot "Lattice strain Q vs pressure" "HydroStrainVsP"
}

##################################################
# plotEpsHStressvsEps33
#	Object: creates a plot of "strains in high stress direction" values for each plane vs. eps33
#	intput:
#	output:
# History
#	Created 15/11/2007 S. Merkel
# Comment:
##################################################
proc plotEpsHStressvsEps33 {} {
	global strain
	global nsteps
	global orientations
	global diff
	global nhkl
	global hkl
	set output [open "EpsHStressVsEps33.tmp" "w"]
	for {set i 0} {$i < $nsteps} {incr i} {
		set line "$strain($i,33)"
		for {set j 0} {$j < $nhkl} {incr j} {
			set eps0 $diff($i,$hkl($j),0)
			set eps90 $diff($i,$hkl($j),90)
			set epsHStress $eps0]
			set line "$line\t$epsHStress"
		}
		puts $output $line
	}
	close $output
	set gnuplot "set ylabel \"Maximum strain(hkl)\"\nset xlabel \"Eps33\"\nset key right top"
	for {set j 0} {$j < $nhkl} {incr j} {
		set col [expr $j+2]
		set nsys [expr $j+1]
		if {$nsys == 1} {
			set gnuplot "$gnuplot\nplot \"EpsHStressVsEps33.tmp\" using 1:$col title \"$hkl($j)\""
		} else {
			set gnuplot "$gnuplot, \"EpsHStressVsEps33.tmp\" using 1:$col title \"$hkl($j)\""
		}
	}
	startGnuplot $gnuplot "Maximum strain vs Eps33" "MaxStrainVsEps33"
}

##################################################
# plotEpsHStressvsP
#	Object: creates a plot of "strains in high stress direction"  values for each plane  vs. pressure
#	intput:
#	output:
# History
#	Created 15/11/2007 S. Merkel
# Comment:
##################################################
proc plotEpsHStressvsP {} {
	global stress
	global nsteps
	global orientations
	global diff
	global nhkl
	global hkl
	set output [open "EpsHStressVsP.tmp" "w"]
	for {set i 0} {$i < $nsteps} {incr i} {
		set line "$stress($i,P)"
		for {set j 0} {$j < $nhkl} {incr j} {
			set eps0 $diff($i,$hkl($j),0)
			set eps90 $diff($i,$hkl($j),90)
			set epsHStress [expr $eps0]
			set line "$line\t$epsHStress"
		}
		puts $output $line
	}
	close $output
	set gnuplot "set ylabel \"Maximum strain (hkl)\"\nset xlabel \"P (GPa)\"\nset key right top"
	for {set j 0} {$j < $nhkl} {incr j} {
		set col [expr $j+2]
		set nsys [expr $j+1]
		if {$nsys == 1} {
			set gnuplot "$gnuplot\nplot \"EpsHStressVsP.tmp\" using 1:$col title \"$hkl($j)\""
		} else {
			set gnuplot "$gnuplot, \"EpsHStressVsP.tmp\" using 1:$col title \"$hkl($j)\""
		}
	}
	startGnuplot $gnuplot "Maximum strain vs pressure" "MaxStrainVsP"
}

##################################################
# plotEpsLStressvsEps33
#	Object: creates a plot of "strains in low stress direction" values for each plane vs. eps33
#	intput:
#	output:
# History
#	Created 15/11/2007 S. Merkel
# Comment:
##################################################
proc plotEpsLStressvsEps33 {} {
	global strain
	global nsteps
	global orientations
	global diff
	global nhkl
	global hkl
	set output [open "EpsLStressVsEps33.tmp" "w"]
	for {set i 0} {$i < $nsteps} {incr i} {
		set line "$strain($i,33)"
		for {set j 0} {$j < $nhkl} {incr j} {
			set eps0 $diff($i,$hkl($j),0)
			set eps90 $diff($i,$hkl($j),90)
			set epsLStress $eps90]
			set line "$line\t$epsLStress"
		}
		puts $output $line
	}
	close $output
	set gnuplot "set ylabel \"Minimum strain(hkl)\"\nset xlabel \"Eps33\"\nset key right top"
	for {set j 0} {$j < $nhkl} {incr j} {
		set col [expr $j+2]
		set nsys [expr $j+1]
		if {$nsys == 1} {
			set gnuplot "$gnuplot\nplot \"EpsLStressVsEps33.tmp\" using 1:$col title \"$hkl($j)\""
		} else {
			set gnuplot "$gnuplot, \"EpsLStressVsEps33.tmp\" using 1:$col title \"$hkl($j)\""
		}
	}
	startGnuplot $gnuplot "Minimum strain vs Eps33" "MinStrainVsEps33"
}

##################################################
# plotEpsLStressvsP
#	Object: creates a plot of "strains in low stress direction"  values for each plane  vs. pressure
#	intput:
#	output:
# History
#	Created 15/11/2007 S. Merkel
# Comment:
##################################################
proc plotEpsLStressvpP {} {
	global stress
	global nsteps
	global orientations
	global diff
	global nhkl
	global hkl
	set output [open "EpsLStressVsP.tmp" "w"]
	for {set i 0} {$i < $nsteps} {incr i} {
		set line "$stress($i,P)"
		for {set j 0} {$j < $nhkl} {incr j} {
			set eps0 $diff($i,$hkl($j),0)
			set eps90 $diff($i,$hkl($j),90)
			set epsLStress [expr $eps90]
			set line "$line\t$epsLStress"
		}
		puts $output $line
	}
	close $output
	set gnuplot "set ylabel \"Minimum strain (hkl)\"\nset xlabel \"P (GPa)\"\nset key right top"
	for {set j 0} {$j < $nhkl} {incr j} {
		set col [expr $j+2]
		set nsys [expr $j+1]
		if {$nsys == 1} {
			set gnuplot "$gnuplot\nplot \"EpsLStressVsP.tmp\" using 1:$col title \"$hkl($j)\""
		} else {
			set gnuplot "$gnuplot, \"EpsLStressVsP.tmp\" using 1:$col title \"$hkl($j)\""
		}
	}
	startGnuplot $gnuplot "Minimum strain vs pressure" "MinStrainVsP"
}

###################################################################################################
#
# Plot strains: open a second window to select the step number and hkl to plot lattice strains for
#
###################################################################################################

##################################################
# plotLatticeStrainUI
#	Object: creates a window to select which reflection and which step number to plot 
#        d-spacings for
# History
#	Created 2/10/2007 S. Merkel
# Comment:
##################################################

proc plotLatticeStrain {widgetindex plotwhat} {
	global nhkl
	global hkl
	global angles
	global nangles
	global diff
	global plotAorCos
	# plot vs angle or 1-3cos2 angle given by $plotAorCos($plotwhat)
	#  if set to cos, plot vis 1-3cos2 psi
	#  else plot vs. angle

	# Plot what?
	set plotsteps [$widgetindex.listboxes.st.steps curselection]
	set plothkl [$widgetindex.listboxes.hkl.hkl curselection]
	# puts "Plotting $plothkl for $plotsteps"
	set nhklplot [llength $plothkl]
	set nstepplot [llength $plotsteps]
	if {($nhklplot == 0) || ($nstepplot == 0)} {return}

	# preparing data to plot
	set str "#"
	for {set i 0} {$i < $nhklplot} {incr i} {
		set thishkl [lindex $plothkl $i]
		set hklname $hkl($thishkl)
		for {set j 0} {$j < $nstepplot} {incr j} {
			set thisstep [lindex $plotsteps $j]
			set str "$str\t$hklname-$thisstep"
		}
	}
	for {set k 0} {$k < $nangles} {incr k} {
		set thisangle $angles($k)
		if {$plotAorCos($plotwhat) == "cos"} {
			set cos [expr cos($thisangle*acos(-1.)/180.)]
			set x [expr 1.-3.*$cos*$cos]
		} else {
			set x $thisangle
		}
		set str "$str\n$x"
		for {set i 0} {$i < $nhklplot} {incr i} {
			set thishkl [lindex $plothkl $i]
			set hklname $hkl($thishkl)
			for {set j 0} {$j < $nstepplot} {incr j} {
				set thisstep [lindex $plotsteps $j]
				set str "$str\t$diff($thisstep,$hklname,$thisangle)"
			}
		}
	}
	# Creating plotting data file
	set output [open "latticestrain.tmp" "w"]
	puts $output $str
	close $output

	# Plotting control
	if {$plotAorCos($plotwhat) == "cos"} {
		set gnuplot "set ylabel \"Lattice strain\"\nset xlabel \"1-3cos^2(angle)\"\nset key right bottom"
		set gnuplot "$gnuplot\nf(x) = a + b * x"
	} else {
		set gnuplot "set ylabel \"Lattice strain\"\nset xlabel \"angle (degrees)\"\nset key right bottom"
		set gnuplot "$gnuplot\nf(x) = a * (1. + b * (1.-3.*(cos(x*pi/180.))**2))"
	}
	set gnuplot "$gnuplot\na = -0.5"
	set gnuplot "$gnuplot\nb = 0.0001"
	set col 1
	for {set i 0} {$i < $nhklplot} {incr i} {
		set thishkl [lindex $plothkl $i]
		for {set j 0} {$j < $nstepplot} {incr j} {
			set thisstep [lindex $plotsteps $j]
			set hklname $hkl($thishkl)
			set thislabel "$hklname$thisstep"
			incr col
			# Fitting removed (errors between GP and TCL)
			set gnuplot "$gnuplot\nfit f(x) \"latticestrain.tmp\" using 1:$col via a, b"
			set gnuplot "$gnuplot\na$thislabel = a"
			set gnuplot "$gnuplot\nb$thislabel = b"
		}
	}
	set col 1
	for {set i 0} {$i < $nhklplot} {incr i} {
		set thishkl [lindex $plothkl $i]
		set hklname $hkl($thishkl)
		for {set j 0} {$j < $nstepplot} {incr j} {
			set thisstep [lindex $plotsteps $j]
			set str "$hklname-$thisstep"
			set thislabel "$hklname$thisstep"
			incr col
			if {$col == 2} {
				set gnuplot "$gnuplot\nplot \"latticestrain.tmp\" using 1:$col with points lt $col ps 1.5 title \"$str\""
				set gnuplot "$gnuplot, a=a$thislabel, b=b$thislabel, f(x) notitle with lines lt $col"
			} else {
				set gnuplot "$gnuplot, \"latticestrain.tmp\" using 1:$col  with points lt $col  ps 1.5 title \"$str\""
				set gnuplot "$gnuplot, a=a$thislabel, b=b$thislabel, f(x) notitle with lines lt $col"
			}
		}
	}
	set gnuplot "$gnuplot;"
	set output [open "plotlatticestrain.tmp" "w"]
	puts $output $gnuplot
	close $output
	startGnuplot $gnuplot "Lattice strains" "latticestrain"

}

proc plotLatticeStrainUI {} {
	global stress
	global strain
	global nsteps
	global nhkl
	global hkl
	
	global indexShowData
	incr indexShowData

	# this variable is passed to know what to plot...
	global plotAorCos
	
	set windowtitle "Plot lattice strain window"
	toplevel .s${indexShowData}
    wm title .s${indexShowData} $windowtitle
	frame .s${indexShowData}.listboxes -bd 1 -relief solid
	# List of available steps
	frame .s${indexShowData}.listboxes.st
	listbox .s${indexShowData}.listboxes.st.steps -height 20 -yscrollcommand ".s${indexShowData}.listboxes.st.stepsscroll set" -selectmode extended -background white -exportselection false
	for {set i 0} {$i < $nsteps} {incr i} {
		set pressure  [format "%4.2f" $stress($i,P)]
		set strain33  [format "%4.3f" $strain($i,33)]
		.s${indexShowData}.listboxes.st.steps insert end "step :$i, $pressure GPa, $strain33"
	}
	scrollbar .s${indexShowData}.listboxes.st.stepsscroll -command ".s${indexShowData}.listboxes.st.steps yview"
	pack .s${indexShowData}.listboxes.st.stepsscroll -side right -fill y
	pack .s${indexShowData}.listboxes.st.steps -padx 5 -pady 5 -side left
	# List of available reflections
	frame .s${indexShowData}.listboxes.hkl
	listbox .s${indexShowData}.listboxes.hkl.hkl -height 20 -yscrollcommand ".s${indexShowData}.listboxes.hkl.hklscroll set" -selectmode extended -background white -exportselection false
	for {set j 0} {$j < $nhkl} {incr j} {
		.s${indexShowData}.listboxes.hkl.hkl insert end "$hkl($j)"
	}
	scrollbar .s${indexShowData}.listboxes.hkl.hklscroll -command ".s${indexShowData}.listboxes.hkl.hkl yview"
	pack .s${indexShowData}.listboxes.hkl.hklscroll -side right -fill y
	pack .s${indexShowData}.listboxes.hkl.hkl -padx 5 -pady 5 -side left
	# Packing lists
	pack .s${indexShowData}.listboxes.st -side left
	pack .s${indexShowData}.listboxes.hkl -side left
	pack .s${indexShowData}.listboxes -pady 5 
	# Angle options
	frame .s${indexShowData}.angles -bd 1 -relief solid
	label .s${indexShowData}.angles.lab -text "Plot vs."
	set plotAorCos(s${indexShowData}plotwhat) "angle"
	radiobutton .s${indexShowData}.angles.angle -text "psi" -anchor w -variable plotAorCos(s${indexShowData}plotwhat) -value angle
	radiobutton .s${indexShowData}.angles.cos -text "(1-3 cos^2(psi))" -anchor w -variable plotAorCos(s${indexShowData}plotwhat) -value cos
	pack .s${indexShowData}.angles.lab .s${indexShowData}.angles.angle .s${indexShowData}.angles.cos -side left
	pack .s${indexShowData}.angles  -pady 5 -fill x
	# Buttons
	frame .s${indexShowData}.buttons
	button .s${indexShowData}.buttons.close -text "Close"  -width 15 \
			-command "destroy  \".s${indexShowData}\""
	button .s${indexShowData}.buttons.export -text "Plot"  -width 15 \
			-command "plotLatticeStrain  \".s${indexShowData}\" \"s${indexShowData}plotwhat\""
    pack .s${indexShowData}.buttons.export -padx 5 -pady 5 -side left
    pack .s${indexShowData}.buttons.close -padx 5 -pady 5 -side left
    pack .s${indexShowData}.buttons
}

###################################################################################################
#
# Plot compression curve: need a dialog to choose which HKL to plot
# Changed 04/2008 to plot an actual fit of eps hydrostatic
#
###################################################################################################

proc plotCompressionCurve {widgetindex} {
	global nhkl
	global hkl
	global angles
	global nangles
	global diff
	global stress
	global nsteps

	# Plot what?
	set plothkl [$widgetindex.listboxes.hkl.hkl curselection]
	# puts "Plotting $plothkl for $plotsteps"
	set nhklplot [llength $plothkl]
	if {($nhklplot == 0)} {return}

	# preparing data to plot
	set str "#"
	for {set i 0} {$i < $nhklplot} {incr i} {
		set thishkl [lindex $plothkl $i]
		set hklname $hkl($thishkl)
		set str "$str\t$hklname-0\t$hklname-hydro\t$hklname-90"
	}
	for {set i 0} {$i < $nsteps} {incr i} {
		set line "$stress($i,P)"
		for {set j 0} {$j < $nhklplot} {incr j} {
			set thishkl [lindex $plothkl $j]
			set eps0 $diff($i,$hkl($thishkl),0)
			set eps90 $diff($i,$hkl($thishkl),90)
			set psids {}
			for {set k 0} {$k < $nangles} {incr k} {
				set psi $angles($k)
				set d $diff($i,$hkl($j),$angles($k))
				lappend psids [list $psi $d]
			}
			foreach {dp Q} [dpQ $psids] break
			set line "$line\t$eps0\t$dp\t$eps90"
		}
		set str "$str\n$line"
	}
	# Creating plotting data file
	set output [open "compression.tmp" "w"]
	puts $output $str
	close $output

	# Plotting control
	set gnuplot "set ylabel \"Strain\"\nset xlabel \"Pressure\"\nset key right top"
	set col 1
	for {set i 0} {$i < $nhklplot} {incr i} {
		set thishkl [lindex $plothkl $i]
		set hklname $hkl($thishkl)
		set str "$hklname"
		set linetype [expr $i+1]
		incr col
		if {$col == 2} {
			set gnuplot "$gnuplot\nplot \"compression.tmp\" using 1:$col notitle with dots lt $linetype"
			incr col
			set gnuplot "$gnuplot, \"compression.tmp\" using 1:$col  with lines lt $linetype title \"$str\""
			incr col
			set gnuplot "$gnuplot, \"compression.tmp\" using 1:$col notitle with dots lt $linetype"
		} else {
			set gnuplot "$gnuplot, \"compression.tmp\" using 1:$col notitle with dots lt $linetype"
			incr col
			set gnuplot "$gnuplot, \"compression.tmp\" using 1:$col  with lines lt $linetype title \"$str\""
			incr col
			set gnuplot "$gnuplot, \"compression.tmp\" using 1:$col notitle with dots lt $linetype"
		}
	}
	set gnuplot "$gnuplot;"
	startGnuplot $gnuplot "Compression curve" "compression"
}

proc plotCompressionCurveUI {} {
	global nhkl
	global hkl
	
	global indexShowData
	incr indexShowData
	
	set windowtitle "Plot compression curve window"
	toplevel .s${indexShowData}
    wm title .s${indexShowData} $windowtitle
	frame .s${indexShowData}.listboxes -bd 1 -relief solid
	# List of available reflections
	frame .s${indexShowData}.listboxes.hkl
	listbox .s${indexShowData}.listboxes.hkl.hkl -height 20 -yscrollcommand ".s${indexShowData}.listboxes.hkl.hklscroll set" -selectmode extended -background white -exportselection false
	for {set j 0} {$j < $nhkl} {incr j} {
		.s${indexShowData}.listboxes.hkl.hkl insert end "$hkl($j)"
	}
	scrollbar .s${indexShowData}.listboxes.hkl.hklscroll -command ".s${indexShowData}.listboxes.hkl.hkl yview"
	pack .s${indexShowData}.listboxes.hkl.hklscroll -side right -fill y
	pack .s${indexShowData}.listboxes.hkl.hkl -padx 5 -pady 5 -side left
	# Packing lists
	pack .s${indexShowData}.listboxes.hkl -side left
	pack .s${indexShowData}.listboxes -pady 5 
	# Buttons
	frame .s${indexShowData}.buttons
	button .s${indexShowData}.buttons.close -text "Close"  -width 15 \
			-command "destroy  \".s${indexShowData}\""
	button .s${indexShowData}.buttons.export -text "Plot"  -width 15 \
			-command "plotCompressionCurve  \".s${indexShowData}\""
    pack .s${indexShowData}.buttons.export -padx 5 -pady 5 -side left
    pack .s${indexShowData}.buttons.close -padx 5 -pady 5 -side left
    pack .s${indexShowData}.buttons
}

###################################################################################################
#
# Display data procedudes: procedures used to display numerical results of EPSC calculations
#
###################################################################################################

##################################################
# exportTxt
#	Object: export data from a showData window to a file
#	intput: 
#		- show data ID
#	output:
# History
#	Created 10/11/2006 S. Merkel
# Comment:
##################################################
proc exportTxt {showWidget} {
	set txt [$showWidget.txt.txt get 1.0 end]
	set file_types {
		{ "Data files" {.dat .DAT} }
		{ "Text files" {.txt .TXT} }
		{ "All files" * }
	}
	set filename [tk_getSaveFile -filetypes $file_types -title "Select file name..." -parent $showWidget]
	if {$filename != ""} {
		set output [open $filename "w"]
		puts $output $txt
		close $output
	}
}

##################################################
# showData
#	Object: creates a txt window to show actual datas
#	intput:
#		- intro: small text describing the data
#		- ncols
#		- nrows
#		- legends: list with a legend for each column
#		- thisdata: array of data (already formatted)
#	output:
# History
#	Created 10/11/2006 S. Merkel
# Comment:
##################################################
proc showData {intro ncols nrows legends thisdata} {
	upvar $thisdata data
	global indexShowData
	incr indexShowData
	set windowtitle [join $legends " "]
	toplevel .s${indexShowData}
    wm title .s${indexShowData} $windowtitle
	frame .s${indexShowData}.txt -bd 4 -relief sunken
	text .s${indexShowData}.txt.txt -yscrollcommand ".s${indexShowData}.txt.scroll set" -background "#FFFFFF" -width 50 -height 30
	scrollbar .s${indexShowData}.txt.scroll -command ".s${indexShowData}.txt.txt yview"
	pack .s${indexShowData}.txt.scroll -side right -fill y
	pack .s${indexShowData}.txt.txt -side left -fill both -expand 1
	pack .s${indexShowData}.txt -fill both -expand 1
	frame .s${indexShowData}.buttons
	button .s${indexShowData}.buttons.close -text "Close"  -width 15 \
			-command "destroy  \".s${indexShowData}\""
	button .s${indexShowData}.buttons.export -text "Export"  -width 15 \
			-command "exportTxt  \".s${indexShowData}\""
    pack .s${indexShowData}.buttons.close -padx 5 -pady 5 -side left
    pack .s${indexShowData}.buttons.export -padx 5 -pady 5 -side left
    pack .s${indexShowData}.buttons
	set txt "# $intro"
	for {set j 0} {$j < $ncols} {incr j} {
		set tmptxt [lindex $legends $j]
		if {$j == 0} {
			set txt "$txt\n#$tmptxt"
		} else {
			set txt "$txt\t$tmptxt"
		}
	}
	set txt "$txt\n"
	for {set i 0} {$i < $nrows} {incr i} {
		for {set j 0} {$j < $ncols} {incr j} {
			if {$j == 0} {
				set txt "$txt$data($i,$j)"
			} else {
				set txt "$txt\t$data($i,$j)"
			}
		}
		set txt "$txt\n"
	}
	.s${indexShowData}.txt.txt insert 1.0 $txt
}


##################################################
# showPTEps33vsStep
#	Object: show the values of eps33, pressure and differential stresses as a funtion of step number
#	intput:
#	output:
# History
#	Created 10/11/2006 S. Merkel
# Comment:
##################################################
proc showPTEps33vsStep {} {
	global nsteps
	global stress
	global strain
	for {set i 0} {$i < $nsteps} {incr i} {
		set data($i,0) $i
		set data($i,1) [format "%1.4f" $strain($i,33)]
		set data($i,2) [format "%3.2f" $stress($i,P)]
		set data($i,3) [format "%1.2f"  $stress($i,t)]
	}
	set legend [list "step" "eps33" "P (GPa)" "t (GPa)"]
	set ncols 4
	set desc "Step, vertical strain, pressure, and differential stresses"
	showData $desc 4 $nsteps $legend data
}

##################################################
# showPEps33ActvsStep
#	Object: show the values of eps33, pressure and slip system activities as a funtion of step number
#	intput:
#	output:
# History
#	Created 10/11/2006 S. Merkel
# Comment:
##################################################
proc showPEps33ActvsStep {} {
	global nsteps
	global stress
	global strain
	global activities
	global nsystems
	for {set i 0} {$i < $nsteps} {incr i} {
		set data($i,0) $i
		set data($i,1) [format "%1.4f" $strain($i,33)]
		set data($i,2) [format "%3.2f" $stress($i,P)]
		for {set j 0} {$j < $nsystems} {incr j} {
			set col [expr $j+3]
			set data($i,$col) [format "%1.3f" $activities($i,$j)]
		}
	}
	set legend [list "step" "eps33" "P (GPa)"]
	for {set j 0} {$j < $nsystems} {incr j} {
		set n [expr $j+1]
		set legend [concat $legend "Sys$n"]
	}
	set ncols [expr 3+$nsystems]
	set desc "Step, vertical strain, pressure, slip system activities"
	showData $desc $ncols $nsteps $legend data
}


###################################################################################################
#
# Main procedudes: load data, build UI...
#
###################################################################################################

##################################################
# loaddata
#	Object: loads data from EVPSC output files
#	intput: difffile: diffraction file
#			phase phase number
#	output:
# History
#	Created 10/10/2006 S. Merkel
#   Changed 07/11/2014 S. Merkel to adapt to EVPSC
# Comment:
##################################################
proc loaddata {difffile phase} {
	countsteps
	readstresses
	readactivities $phase
	readdiffraction $difffile $phase
}

##################################################
# clearTmpFiles
#	Object: delete temporary files created for plotting stuff
#	intput:
#	output:
# History
#	2016-10: added QvsStep
#	Created 10/11/2006 S. Merkel
# Comment:
##################################################
proc clearTmpFiles {} {
	file delete "pVsEps33.tmp"
	file delete "tVsEps33.tmp"
	file delete "QVsEps33.tmp"
	file delete "stressVsEps33.tmp"
	file delete "SigmaQVsEps33.tmp"
	file delete "EpsHydroVsEps33.tmp"
	file delete "EpsHStressVsEps33.tmp"
	file delete "EpsLStressVsEps33.tmp"
	file delete "stressVsEps33.tmp"
	file delete "tVsP.tmp"
	file delete "QVsP.tmp"
	file delete "QVsStep.tmp"
	file delete "SigmaQVsP.tmp"
	file delete "EpsHydroVsP.tmp"
	file delete "EpsHStressVsP.tmp"
	file delete "EpsLStressVsP.tmp"
	file delete "QVsT.tmp"
	file delete "actVsEps33.tmp"
	file delete "actVsStep.tmp"
	file delete "actVsP.tmp"
	file delete "tmp.eps"
	file delete "latticestrain.tmp"
	file delete "resultat.tk"
	file delete "compression.tmp"
	file delete "plotlatticestrain.tmp"
	file delete "fit.log"
}

##################################################
# clearImages
#	Object: delete image files create previously (EPS and GIF)
#	intput:
#	output:
# History
#	2016-10: added QvsStep
#	Created 10/11/2006 S. Merkel
# Comment:
##################################################
proc clearImages {} {
	file delete "pVsEps33.gif"
	file delete "QVsEps33.gif"
	file delete "tVsEps33.gif"
	file delete "tVsP.gif"
	file delete "QVsP.gif"
	file delete "QVsT.gif"
	file delete "QVsStep.gif"
	file delete "actVsEps33.gif"
	file delete "actVsStep.gif"
	file delete "actVsP.gif"
	file delete "latticestrain.gif"
	file delete "pVsEps33.eps"
	file delete "QVsEps33.eps"
	file delete "tVsEps33.eps"
	file delete "tVsP.eps"
	file delete "QVsP.eps"
	file delete "QVsT.eps"
	file delete "QVsStep.eps"
	file delete "actVsEps33.eps"
	file delete "actVsStep.eps"
	file delete "actVsP.eps"
	file delete "latticestrain.eps"
}

##################################################
# buildUI
#	Object: builds the user interface
#	intput:
#	output:
# History
#	2016-10: added QvsStep
#	Created 10/10/2006 S. Merkel
# Comment:
##################################################
proc buildUI {} {
	# Set title for main window
	wm title . "Analysis or EPSC results"
	# Plot results frame
	frame .buttonplots -borderwidth 1 -relief solid
	frame .buttonplots.b
	label .buttonplots.b.l1 -text "Plots"
	button .buttonplots.b.plotPvsEps -text "P vs Eps33"  -width 20 -command {plotPvsEps33}
	button .buttonplots.b.plotTvsEps -text "t vs Eps33"  -width 20 -command {plotTvsEps33}
	button .buttonplots.b.plotStressvsEps -text "Sigma vs Eps33"  -width 20 -command {plotStressesvsEps33}
	button .buttonplots.b.plotAvsEps -text "Act. vs Eps33"  -width 20 -command {plotAvsEps33}
	button .buttonplots.b.plotQvsEps -text "Q vs Eps33"  -width 20 -command {plotQvsEps33}
	button .buttonplots.b.plotErrQvsEps -text "Sigma(Q) vs Eps33"  -width 20 -command {plotErrQvsEps33}
	button .buttonplots.b.plotEpsHydrovsEps -text "Hydro. strain vs Eps33"  -width 20 -command {plotEpsHydrovsEps33}
	button .buttonplots.b.plotEpsHStressvsEps -text "Max. strain vs Eps33"  -width 20 -command {plotEpsHStressvsEps33}
	button .buttonplots.b.plotEpsLStressvsEps -text "Min. strain vs Eps33"  -width 20 -command {plotEpsLStressvsEps33}
	button .buttonplots.b.plotTvsP -text "t vs P"  -width 20 -command {plotTvsP}
	button .buttonplots.b.plotAvsP -text "Act. vs P"  -width 20 -command {plotAvsP}
	button .buttonplots.b.plotQvsP -text "Q vs P"  -width 20 -command {plotQvsP}
	button .buttonplots.b.plotErrQvsP -text "Sigma(Q) vs P"  -width 20 -command {plotErrQvsP}
	button .buttonplots.b.plotEpsHydrovsP -text "Hydro. strain vs P"  -width 20 -command {plotEpsHydrovsP}
	button .buttonplots.b.plotEpsHStressvsP -text "Max. strain vs P"  -width 20 -command {plotEpsHStressvsP}
	button .buttonplots.b.plotEpsLStressvsP -text "Min. strain vs P"  -width 20 -command {plotEpsLStressvpP}
	button .buttonplots.b.plotAvsStep -text "Act. vs step"  -width 20 -command {plotAvsStep}
	button .buttonplots.b.plotQvsStep -text "Q vs step"  -width 20 -command {plotQvsStep}
	button .buttonplots.b.plotQvsT -text "Q vs t"  -width 20 -command {plotQvsT}
	button .buttonplots.b.latticestrain -text "Lattice strain"  -width 20 -command {plotLatticeStrainUI}
	button .buttonplots.b.compressioncurve -text "Compression Curve"  -width 20 -command {plotCompressionCurveUI}
	label .buttonplots.b.l2 -text "Data"
	button .buttonplots.b.showPTEps33vsStep -text "P, t, Eps33 vs Step"  -width 20 -command {showPTEps33vsStep}
	button .buttonplots.b.showPEps33ActvsStep -text "P, Eps33, Act vs Step"  -width 20 -command {showPEps33ActvsStep}

	grid config .buttonplots.b.l1 -column 1 -row 0 -padx 5 -pady 5
	grid config .buttonplots.b.l2 -column 3 -row 0 -padx 5 -pady 5

	grid config .buttonplots.b.plotPvsEps -column 0 -row 1 -padx 5 -pady 5
	grid config .buttonplots.b.plotTvsEps -column 0 -row 2 -padx 5 -pady 5
	grid config .buttonplots.b.plotStressvsEps -column 0 -row 3 -padx 5 -pady 5
	grid config .buttonplots.b.plotAvsEps -column 0 -row 4 -padx 5 -pady 5
	grid config .buttonplots.b.plotQvsEps -column 0 -row 5 -padx 5 -pady 5
	grid config .buttonplots.b.plotErrQvsEps -column 0 -row 6 -padx 5 -pady 5
	grid config .buttonplots.b.plotEpsHydrovsEps -column 0 -row 7 -padx 5 -pady 5
	grid config .buttonplots.b.plotEpsHStressvsEps -column 0 -row 8 -padx 5 -pady 5
	grid config .buttonplots.b.plotEpsLStressvsEps -column 0 -row 9 -padx 5 -pady 5
	grid config .buttonplots.b.plotTvsP -column 1 -row 1 -padx 5 -pady 5
	grid config .buttonplots.b.plotAvsP -column 1 -row 2 -padx 5 -pady 5
	grid config .buttonplots.b.plotQvsP -column 1 -row 3 -padx 5 -pady 5
	grid config .buttonplots.b.plotErrQvsP -column 1 -row 4 -padx 5 -pady 5
	grid config .buttonplots.b.plotEpsHydrovsP -column 1 -row 5 -padx 5 -pady 5
	grid config .buttonplots.b.plotEpsHStressvsP -column 1 -row 6 -padx 5 -pady 5
	grid config .buttonplots.b.plotEpsLStressvsP -column 1 -row 7 -padx 5 -pady 5
	grid config .buttonplots.b.plotAvsStep  -column 2 -row 1 -padx 5 -pady 5
	grid config .buttonplots.b.plotQvsStep -column 2 -row 2 -padx 5 -pady 5
	grid config .buttonplots.b.plotQvsT -column 2 -row 3 -padx 5 -pady 5
	grid config .buttonplots.b.latticestrain -column 2 -row 3 -padx 5 -pady 5
	grid config .buttonplots.b.compressioncurve -column 2 -row 4 -padx 5 -pady 5

	grid config .buttonplots.b.showPTEps33vsStep -column 3 -row 1 -padx 5 -pady 5
	grid config .buttonplots.b.showPEps33ActvsStep -column 3 -row 2 -padx 5 -pady 5

	grid config .buttonplots.b
	# Erase files frame
	frame .buttonsErase -borderwidth 1 -relief solid
	label .buttonsErase.l -text "Clear"
	frame .buttonsErase.b
	button .buttonsErase.b.cleartmp -text "Clear tmp files"  -width 20 -command {clearTmpFiles}
	button .buttonsErase.b.clearimages -text "Clear images"  -width 20 -command {clearImages}
    pack .buttonsErase.b.cleartmp -padx 5 -pady 5 -side left
	pack .buttonsErase.b.clearimages -padx 5 -pady 5 -side left
	pack .buttonsErase.l
	pack  .buttonsErase.b
	# End application
	label .titre -text " Analysis of results from EVPSC calculations " -bd 2 -relief ridge -padx 5 -pady 5
	button .done -text "Exit"  -width 10 -command {exit}
	label .copyright -text "21 feb 2017, version 3.1 - 2006-2017, S. Merkel, Universite Lille, France"
	# Finishing up
    pack .titre -padx 5 -pady 5
    pack .buttonplots -padx 5 -pady 5
	pack .buttonsErase -padx 5 -pady 5
    pack .done -padx 5 -pady 5
	pack .copyright
}


###################################################################################################
#
# Main
#
###################################################################################################

if { $argc != 2 } {
        puts "This script requires exactly 2 arguments. The first argument is the name of the diffraction file defined in evpsc.in. The second argument is the phase number"
        puts "For example, ./parseEVPSC.tcl fcc.dif2 1"
        puts "Please try again."
		exit
} else {
	set difffile [lindex $argv 0]
	set phase [lindex $argv 1]
	set indexPlot 0
	set indexShowData 0
	loaddata $difffile $phase
	buildUI
}