# Fidelity Pointwise V18.6R3 Journal file - Sun Jan 14 18:26:08 2024

package require PWI_Glyph 6.22.1

pw::Application setUndoMaximumLevels 5
pw::Application reset
pw::Application markUndoLevel {Journal Reset}

pw::Application clearModified

#################################################################
###########################OPTIONS###############################

set nDim 2

set gridType "NS"
set GRBL 1.15
set BLSpacing 8e-6

set Angle 5

set path "/home/pagaf/Desktop/masterThesis/Multielement/Meshes/Mesh_"

set iter 6
set iterm1 [expr {$iter - 1}]
# set iterm1 "0_ModifiedLETE"

set path2Mesh $path$iter
set extension ".pw"
set ReferenceMeshName $path$iterm1$extension

set Boundaries2Ref [list "Slat" "Airfoil" "Flap"]
set Boundaries2Preserve [list "Slat" "Airfoil" "Flap"]

set existingDomainsOrBlocks [list "dom-1"]
set isNSOrEuler [list "NS"]
set MaxEdge [list 4.0]

#################################################################
#################################################################

proc MergeLists { List2TakeFrom List2PutInto } {


  # Check if connector is already present
  for {set j 0} {$j < [llength $List2TakeFrom]} {incr j} {

    set isIn false

    for {set k 0} {$k < [llength $List2PutInto] && !$isIn} {incr k} {
      if {[lindex $List2TakeFrom $j] == [lindex $List2PutInto $k]} {
        set isIn true
      }
    }

    if {!$isIn} {
      lappend List2PutInto [lindex $List2TakeFrom $j]
    }

  }

  return $List2PutInto

}



proc RemoveFromList { List2Remove List2Clean } {

  set ClearedList [list]

  # Check if connector is already present
  for {set j 0} {$j < [llength $List2Clean]} {incr j} {

    set isIn false

    for {set k 0} {$k < [llength $List2Remove] && !$isIn} {incr k} {
      if {[lindex $List2Clean $j] == [lindex $List2Remove $k]} {
        set isIn true
      }
    }

    if {!$isIn} {
      lappend ClearedList [lindex $List2Clean $j]
    }

  }

  return $ClearedList

}

proc RemoveDuplicatesFromList { List2Clean } {

  set ClearedList [list]

  # Check if connector is already present
  for {set j 0} {$j < [llength $List2Clean]} {incr j} {

    set isIn false

    for {set k 0} {$k < [llength $ClearedList] && !$isIn} {incr k} {
      if {[lindex $List2Clean $j] == [lindex $ClearedList $k]} {
        set isIn true
      }
    }

    if {!$isIn} {
      lappend ClearedList [lindex $List2Clean $j]
    }

  }

  return $ClearedList

}


#################################################################

puts "Opening file..."

pw::Application reset -keep Clipboard
set _TMP(mode_1) [pw::Application begin ProjectLoader]
  $_TMP(mode_1) initialize $ReferenceMeshName
  $_TMP(mode_1) setAppendMode false
  $_TMP(mode_1) setRepairMode Defer
  $_TMP(mode_1) load
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application resetUndoLevels

pw::Connector setCalculateDimensionMaximum 100000

set Volumes [list]
set Domains [list]

if { $nDim == 2 } {

  foreach domainName $existingDomainsOrBlocks {
    lappend Domains [pw::GridEntity getByName $domainName]
  }


} else {
  puts "Emptying existing volumes..."
  foreach blockName $existingDomainsOrBlocks {
    lappend Volumes [pw::GridEntity getByName $blockName]
  }
  set _TMP(mode_1) [pw::Application begin UnstructuredSolver $Volumes]
    $_TMP(mode_1) run Release
  $_TMP(mode_1) end
  unset _TMP(mode_1)
  pw::Application markUndoLevel Solve
}


puts "Importing point cloud..."


set cloud [pw::SourcePointCloud create]
set srcPts [list]

set a [open "./PointCloud.dat"]
set lines [split [read $a] "\n"]
close $a;                          # Saves a few bytes :-)
foreach line $lines {

  set wordList [regexp -inline -all -- {\S+} $line]
  # puts $wordList

  if { [llength $wordList] > 1 } {
    set newPoint [pwu::Vector3 set [lindex $wordList 0] [lindex $wordList 1] [lindex $wordList 2]]
    lappend srcPts [list $newPoint [lindex $wordList 3] [lindex $wordList 4]]
  }
}

puts "Saving connectors distributions..."

set Conns2Preserve [list]
set Doms2Ref [list]

if { $nDim == 2 } {

  foreach Boundary $Boundaries2Preserve {

    set BC [pw::BoundaryCondition getByName $Boundary]
    set Cons [$BC getEntities]
    foreach con $Cons {
      lappend Conns2Preserve $con
    }

  }

} else {
  # First of all assure that the current spacings won't be coarsened more
  foreach Boundary $Boundaries2Preserve {

    set BC [pw::BoundaryCondition getByName $Boundary]
    set Domains [$BC getEntities]
    foreach Domain $Domains {
      foreach edge [$Domain getEdges] {
        foreach con [$edge getConnectors] {
          lappend Conns2Preserve $con
        }
      }
    }

  }
}


set Conns2Preserve [RemoveDuplicatesFromList $Conns2Preserve]

foreach con2Preserve $Conns2Preserve {

  # puts  [$con2Preserve getName]

  set NPoints [$con2Preserve getDimensions]
  for {set i 1} {$i < $NPoints} {incr i} {

    set Point [$con2Preserve getXYZ -grid $i]
    set PointAfter [$con2Preserve getXYZ -grid [expr {$i+1}]]

    set Difference [pwu::Vector3 subtract $PointAfter $Point]
    set Distance [pwu::Vector3 length $Difference]
    lappend srcPts [list $Point $Distance 0.85]
  }
}


$cloud addPoints $srcPts


set Conns2Ref [list]
if { $nDim == 2 } {
  set Conns2Ref [list]
  foreach Boundary $Boundaries2Ref {

    set BC [pw::BoundaryCondition getByName $Boundary]
    set Cons [$BC getEntities]
    foreach con $Cons {
      lappend Conns2Ref $con
    }
  }

} else {
  # Then select all of the connectors from the BCs to adapt
  foreach Boundary $Boundaries2Ref {

    set BC [pw::BoundaryCondition getByName $Boundary]
    set Domains [$BC getEntities]
    foreach Domain $Domains {
      foreach edge [$Domain getEdges] {
        foreach con [$edge getConnectors] {
          lappend Conns2Ref $con
        }
      }
    }
  }

}


set Conns2Ref [RemoveDuplicatesFromList $Conns2Ref]


set StructCons [list]
set StructPairingCons [list]
set ParallelCons [list]
set ToConsider [list]

if { $nDim == 3 } {
  # From these I have to remove all of the connectors that belong to structured domains
  set StructDoms [pw::Grid getAll -type pw::DomainStructured]
  # Now cycle on all of the struct domains and create a list of all of their connectors

  foreach structDom $StructDoms {

    set Edges [$structDom getEdges]

    foreach Edge $Edges {
      set Cons [$Edge getConnectors]
      foreach Con $Cons {
        lappend StructCons $Con
      }
    }


    # There are 4 edges
    set Edge [lindex $Edges 0]
    set Cons [$Edge getConnectors]
    set Edge [lindex $Edges 2]
    set Cons2 [$Edge getConnectors]
    lappend StructPairingCons [list $Cons $Cons2]

    set Edge [lindex $Edges 1]
    set Cons [$Edge getConnectors]
    set Edge [lindex $Edges 3]
    set Cons2 [$Edge getConnectors]
    lappend StructPairingCons [list $Cons $Cons2]


  }

  # Now I have to rearrange all of the connectors such that I have all of the parallel connectors
  # This is all based on the fact that only one connector is present per each edge
  set ParallelCons [list]
  set ToConsider [list]
  for {set i 0} {$i < [llength $StructPairingCons]} {incr i} {
    lappend ToConsider 1
  }
  for {set i 0} {$i < [llength $StructPairingCons]} {incr i} {

    if { [lindex $ToConsider $i ] } {
      set pair [lindex $StructPairingCons $i]

      set list2Append $pair

      for {set j [expr {$i+1}]} {$j < [llength $StructPairingCons]} {incr j} {

        if { [lindex $ToConsider $j ] } {

          set pairTwo [lindex $StructPairingCons $j]

          if { [lindex $pair 0] == [lindex $pairTwo 0] } {
            lappend list2Append [lindex $pairTwo 1]
            lset ToConsider $j 0
          }

          if { [lindex $pair 0] == [lindex $pairTwo 1] } {
            lappend list2Append [lindex $pairTwo 0]
            lset ToConsider $j 0
          }

          if { [lindex $pair 1] == [lindex $pairTwo 0] } {
            lappend list2Append [lindex $pairTwo 1]
            lset ToConsider $j 0
          }

          if { [lindex $pair 1] == [lindex $pairTwo 1] } {
            lappend list2Append [lindex $pairTwo 0]
            lset ToConsider $j 0
          }

        }

      }

      lappend ParallelCons $list2Append

    }
  }

  set StructCons [RemoveDuplicatesFromList $StructCons]

}

puts "Refining standard connectors... "

if { $nDim == 2 } {
  pw::DomainUnstructured setInitializeInterior 0
}

foreach Con2Ref $Conns2Ref {

  set _TMP(mode_1) [pw::Application begin Modify [list $Con2Ref]]
    pw::Connector swapDistribution Tanh [list [list $Con2Ref 1]]
  $_TMP(mode_1) end
  unset _TMP(mode_1)
  pw::Application markUndoLevel Distribute

}

# Now get all of the clouds
set clouds2Use [pw::Source getAll]
pw::Connector setDimensionFromSizeField -include $clouds2Use -calculationMethod MinimumValue $Conns2Ref


if { $nDim == 3 } {
  puts "Refining connectors from structured domains... "
  foreach cons $ParallelCons {
    foreach Con2Ref $cons {

      set _TMP(mode_1) [pw::Application begin Modify [list $Con2Ref]]
        pw::Connector swapDistribution Tanh [list [list $Con2Ref 1]]
      $_TMP(mode_1) end
      unset _TMP(mode_1)
      pw::Application markUndoLevel Distribute

    }
    pw::Connector setDimensionFromSizeField -include $clouds2Use -matchDimensions -calculationMethod MinimumValue $cons
  }
}

if { $nDim == 2 } {
  pw::DomainUnstructured setInitializeInterior 1
}

set Doms2Ref [list]
if { $nDim == 3 } {
  puts "Refining selected domains..."
  foreach Boundary $Boundaries2Ref {

    set BC [pw::BoundaryCondition getByName $Boundary]
    set Domains [$BC getEntities]
    foreach Domain $Domains {
      lappend Doms2Ref $Domain
      puts [$Domain getName]
    }
  }

  set Doms2Ref [RemoveFromList $StructDoms $Doms2Ref]

  if { [llength $Doms2Ref] > 0 } {
    set _TMP(mode_1) [pw::Application begin UnstructuredSolver $Doms2Ref]
      $_TMP(mode_1) run Initialize
    $_TMP(mode_1) end
    unset _TMP(mode_1)
    pw::Application markUndoLevel Solve
  }
}


if { $nDim == 2 } {

  puts "Setting Domain quantities..."

  for {set i 0} {$i < [llength $Domains]} {incr i} {
    set Domain [lindex $Domains $i]
    set _TMP(mode_1) [pw::Application begin UnstructuredSolver [list $Domain]]

      if { [lindex $isNSOrEuler $i] == "NS" } {
        $Domain setUnstructuredSolverAttribute TRexGrowthRate $GRBL
        set _TMP(PW_1) [pw::TRexCondition getByName Wall]
        $_TMP(PW_1) setValue $BLSpacing
        unset _TMP(PW_1)
      }

      $Domain setSizeFieldDecay 0.8
      $Domain setUnstructuredSolverAttribute EdgeMaximumLength [lindex $MaxEdge $i]
    $_TMP(mode_1) end
    unset _TMP(mode_1)
    pw::Application markUndoLevel Solve
  }

  for {set i 0} {$i < [llength $Domains]} {incr i} {
    set Domain [lindex $Domains $i]
    set _TMP(mode_1) [pw::Application begin UnstructuredSolver [list $Domain]]

      $_TMP(mode_1) setStopWhenFullLayersNotMet true
      $_TMP(mode_1) setAllowIncomplete true
      $_TMP(mode_1) run Initialize

    $_TMP(mode_1) end
    unset _TMP(mode_1)
    pw::Application markUndoLevel Solve
  }


} else {

  puts "Setting Block quantities..."

  for {set i 0} {$i < [llength $Volumes]} {incr i} {
    set Volume [lindex $Volumes $i]
    set _TMP(mode_1) [pw::Application begin UnstructuredSolver [list $Volume]]

      if { [lindex $isNSOrEuler $i] == "NS" } {
        $Volume setUnstructuredSolverAttribute TRexGrowthRate $GRBL
        set _TMP(PW_1) [pw::TRexCondition getByName Wall]
        $_TMP(PW_1) setValue $BLSpacing
        unset _TMP(PW_1)
      }

      $Volume setSizeFieldDecay 0.8
      $Volume setUnstructuredSolverAttribute EdgeMaximumLength [lindex $MaxEdge $i]
    $_TMP(mode_1) end
    unset _TMP(mode_1)
    pw::Application markUndoLevel Solve
  }

  puts "Initialize blocks..."

  for {set i 0} {$i < [llength $Volumes]} {incr i} {
    set Volume [lindex $Volumes $i]
    set _TMP(mode_1) [pw::Application begin UnstructuredSolver [list $Volume]]
      $_TMP(mode_1) setStopWhenFullLayersNotMet true
      $_TMP(mode_1) setAllowIncomplete true
      $_TMP(mode_1) run Initialize
      $Volume setUnstructuredSolverAttribute WCNSmoothConvergenceCostThreshold 0.9
      $Volume setUnstructuredSolverAttribute WCNSmoothCostAngleThreshold 160
      $Volume setUnstructuredSolverAttribute WCNSmoothRelaxationFactor 0.15
      $_TMP(mode_1) run Smooth 400
    $_TMP(mode_1) end
    unset _TMP(mode_1)
    pw::Application markUndoLevel Initialize
  }

}


set Extension ".su2"

pw::Application setCAESolver SU2 $nDim
pw::Application markUndoLevel {Select Solver}

# set Extension ".pw"
# pw::Application save $path2Mesh$RefLevel$Extension

if { $nDim == 2 } {
  set _TMP(mode_1) [pw::Application begin CaeExport [pw::Entity sort $Domains]]
    $_TMP(mode_1) initialize -strict -type CAE $path2Mesh$Extension
    $_TMP(mode_1) setAttribute FilePrecision Double
    $_TMP(mode_1) verify
    $_TMP(mode_1) write
  $_TMP(mode_1) end
  unset _TMP(mode_1)
} else {
  set _TMP(mode_1) [pw::Application begin CaeExport [pw::Entity sort $Volumes]]
    $_TMP(mode_1) initialize -strict -type CAE $path2Mesh$Extension
    $_TMP(mode_1) setAttribute FilePrecision Double
    $_TMP(mode_1) verify
    $_TMP(mode_1) write
  $_TMP(mode_1) end
  unset _TMP(mode_1)
}

set Extension ".pw"
pw::Application save $path2Mesh$Extension
