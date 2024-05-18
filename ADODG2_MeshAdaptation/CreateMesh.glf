# Fidelity Pointwise V18.6R3 Journal file - Sun Jan 14 18:26:08 2024

package require PWI_Glyph 6.22.1

pw::Application setUndoMaximumLevels 5
pw::Application reset
pw::Application markUndoLevel {Journal Reset}

pw::Application clearModified

source "./CreateMesh_Params.glf"
source "./CreateMesh_Fun.glf"

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

# Removing all of the previous sources
pw::Entity delete [pw::Source getAll]

set cloud [pw::SourcePointCloud create]
set srcPts [list]

set filename "./PointCloud.dat"
if { $BLTreatment } {
  puts "BL Treatment is set to True.. Opening the surface point cloud"
  set filename "./surface_PointCloud.dat"
}

set a [open $filename]
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

set Extension ".pw"
pw::Application save $path2Mesh$Extension


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

# Structured domains avoided for now

set StructDoms [pw::Grid getAll -type pw::DomainStructured]

# set StructCons [list]
# set StructPairingCons [list]
# set ParallelCons [list]
# set ToConsider [list]
#
# if { $nDim == 3 } {
#   # From these I have to remove all of the connectors that belong to structured domains
#   set StructDoms [pw::Grid getAll -type pw::DomainStructured]
#   # Now cycle on all of the struct domains and create a list of all of their connectors
#
#   foreach structDom $StructDoms {
#
#     set Edges [$structDom getEdges]
#
#     foreach Edge $Edges {
#       set Cons [$Edge getConnectors]
#       foreach Con $Cons {
#         lappend StructCons $Con
#       }
#     }
#
#
#     # There are 4 edges
#     set Edge [lindex $Edges 0]
#     set Cons [$Edge getConnectors]
#     set Edge [lindex $Edges 2]
#     set Cons2 [$Edge getConnectors]
#     lappend StructPairingCons [list $Cons $Cons2]
#
#     set Edge [lindex $Edges 1]
#     set Cons [$Edge getConnectors]
#     set Edge [lindex $Edges 3]
#     set Cons2 [$Edge getConnectors]
#     lappend StructPairingCons [list $Cons $Cons2]
#
#
#   }
#
#   # Now I have to rearrange all of the connectors such that I have all of the parallel connectors
#   # This is all based on the fact that only one connector is present per each edge
#   set ParallelCons [list]
#   set ToConsider [list]
#   for {set i 0} {$i < [llength $StructPairingCons]} {incr i} {
#     lappend ToConsider 1
#   }
#   for {set i 0} {$i < [llength $StructPairingCons]} {incr i} {
#
#     if { [lindex $ToConsider $i ] } {
#       set pair [lindex $StructPairingCons $i]
#
#       set list2Append $pair
#
#       for {set j [expr {$i+1}]} {$j < [llength $StructPairingCons]} {incr j} {
#
#         if { [lindex $ToConsider $j ] } {
#
#           set pairTwo [lindex $StructPairingCons $j]
#
#           if { [lindex $pair 0] == [lindex $pairTwo 0] } {
#             lappend list2Append [lindex $pairTwo 1]
#             lset ToConsider $j 0
#           }
#
#           if { [lindex $pair 0] == [lindex $pairTwo 1] } {
#             lappend list2Append [lindex $pairTwo 0]
#             lset ToConsider $j 0
#           }
#
#           if { [lindex $pair 1] == [lindex $pairTwo 0] } {
#             lappend list2Append [lindex $pairTwo 1]
#             lset ToConsider $j 0
#           }
#
#           if { [lindex $pair 1] == [lindex $pairTwo 1] } {
#             lappend list2Append [lindex $pairTwo 0]
#             lset ToConsider $j 0
#           }
#
#         }
#
#       }
#
#       lappend ParallelCons $list2Append
#
#     }
#   }
#
#   set StructCons [RemoveDuplicatesFromList $StructCons]
#
# }

puts "Changing connectors distributions to Tanh... "

pw::DomainUnstructured setInitializeInterior 0

foreach Con2Ref $Conns2Ref {

  set _TMP(mode_1) [pw::Application begin Modify [list $Con2Ref]]
    pw::Connector swapDistribution Tanh [list [list $Con2Ref 1]]
  $_TMP(mode_1) end
  unset _TMP(mode_1)
  pw::Application markUndoLevel Distribute

}

# Refine directly the connectors if 2D
if { $nDim == 2 } {
  # Now get all of the clouds
  pw::Connector setDimensionFromSizeField -include $cloud -calculationMethod MinimumValue $Conns2Ref
}

# Structured domains are avoided for now

# if { $nDim == 3 } {
#   puts "Refining connectors from structured domains... "
#   foreach cons $ParallelCons {
#     foreach Con2Ref $cons {
#
#       set _TMP(mode_1) [pw::Application begin Modify [list $Con2Ref]]
#         pw::Connector swapDistribution Tanh [list [list $Con2Ref 1]]
#       $_TMP(mode_1) end
#       unset _TMP(mode_1)
#       pw::Application markUndoLevel Distribute
#
#     }
#     pw::Connector setDimensionFromSizeField -include $cloud -matchDimensions -calculationMethod MinimumValue $cons
#   }
# }

pw::DomainUnstructured setInitializeInterior 1

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
        set source [pw::SourceEntity getByName src-1]
        foreach dom $Doms2Ref {
          $dom excludeSource $source false
          # Somehow I always need to set this up
          $dom setUnstructuredSolverAttribute TRexPushAttributes True
        }
      $_TMP(mode_1) run Initialize
    $_TMP(mode_1) end
    unset _TMP(mode_1)
    pw::Application markUndoLevel Solve
  }
}


set Extension ".pw"
pw::Application save $path2Mesh$Extension

# Now, if BL treatment was on, I have to remove all of the points from the cloud and create a new one
# with only the points for the volume


if { $BLTreatment } {

  # Removing all of the previous sources
  pw::Entity delete [pw::Source getAll]

  set cloud [pw::SourcePointCloud create]
  set srcPts [list]

  puts "BL Treatment is set to True.. Opening the volume point cloud"
  set filename "./PointCloud.dat"

  set a [open $filename]
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

  $cloud addPoints $srcPts

  # And now adapt only with the volume cloud
}





if { $nDim == 2 } {

  puts "Setting Domain quantities..."

  for {set i 0} {$i < [llength $Domains]} {incr i} {
    set Domain [lindex $Domains $i]
    set _TMP(mode_1) [pw::Application begin UnstructuredSolver [list $Domain]]

    set gridType [lindex $isNSOrEuler $i]

    if { $gridType == "NS" } {

      $Domain setUnstructuredSolverAttribute TRexGrowthRate $GRBL
      set _TMP(PW_1) [pw::TRexCondition getByName $TRexBCName]
      $_TMP(PW_1) setValue $BLSpacing
      unset _TMP(PW_1)

      set minValue 1e5
      set maxValue -1

      foreach NSBound $NSBoundaries {
        set BC [pw::BoundaryCondition getByName $NSBound]
        set Cons [$BC getEntities]
        foreach con $Cons {
          set NPoints [$con getDimensions]
          for {set iPoint 1} {$iPoint < $NPoints} {incr iPoint} {

            set Point [$con getXYZ -grid $iPoint]
            set PointAfter [$con getXYZ -grid [expr {$iPoint+1}]]

            set Difference [pwu::Vector3 subtract $PointAfter $Point]
            set Distance [pwu::Vector3 length $Difference]

            set minValue [expr {min($minValue, $Distance)}]
            set maxValue [expr {max($maxValue, $Distance)}]
          }
        }
      }

      set ToWrite "Max spacing on the TRex BC = "
      puts $ToWrite$maxValue
      set ToWrite "Min spacing on the TRex BC = "
      puts $ToWrite$minValue

      set ToWrite "Number of layers needed for min spacing to reach isotropy = "
      set minNOfLayers [expr {int(log($minValue/$BLSpacing)/log($GRBL))}]
      puts $ToWrite$minNOfLayers

      set ToWrite "Number of layers needed for max spacing to reach set wall distance = "
      set maxNLayers [expr {1- $BLHeight * (1-$GRBL) / $BLSpacing}]
      set maxNLayers [expr {int(log($maxNLayers)/log($GRBL))}]
      puts $ToWrite$maxNLayers

      set NMaxLayers [expr {max($minNOfLayers, $maxNLayers)}]

      $Domain setUnstructuredSolverAttribute TRexMaximumLayers $NMaxLayers
      $Domain setUnstructuredSolverAttribute TRexIsotropicHeight $IsotropicHeight
      $Domain setUnstructuredSolverAttribute TRexFullLayers $NFullLayers

    }

    if { $gridType == "Euler" } {
      $Domain setUnstructuredSolverAttribute TRexMaximumLayers 0
      $Domain setUnstructuredSolverAttribute TRexFullLayers 0
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

    set source [pw::SourceEntity getByName src-1]
    $Domain excludeSource $source false

    set ToWrite "Initializing domain..."
    set ToWrite_1 [lindex $existingDomainsOrBlocks $i]
    set ToWrite_2 "..."
    puts $ToWrite$ToWrite_1$ToWrite_2

    $_TMP(mode_1) run Initialize

    $_TMP(mode_1) end
    unset _TMP(mode_1)
    pw::Application markUndoLevel Solve

    set ToWrite "Domain created. Number of Points = "
    set NPoints [lindex [$Domain getDimensions] 0]
    set ToWrite_2 ". Number of Cells = "
    set NCells [$Domain getCellCount]
    puts $ToWrite$NPoints$ToWrite_2$NCells
  }




} else {

  puts "Setting Block quantities..."

  for {set i 0} {$i < [llength $Volumes]} {incr i} {
    set Volume [lindex $Volumes $i]
    set _TMP(mode_1) [pw::Application begin UnstructuredSolver [list $Volume]]

    set gridType [lindex $isNSOrEuler $i]

    if { $gridType == "NS" } {

      $Volume setUnstructuredSolverAttribute TRexGrowthRate $GRBL
      set _TMP(PW_1) [pw::TRexCondition getByName $TRexBCName]
      $_TMP(PW_1) setValue $BLSpacing
      unset _TMP(PW_1)

      set minValue 1e5
      set maxValue -1

      foreach NSBound $NSBoundaries {
        set BC [pw::BoundaryCondition getByName $NSBound]
        set Doms [$BC getEntities]
        foreach Dom $Doms {
          set NCells [$Dom getCellCount]
          for {set iCell 1} {$iCell <= $NCells} {incr iCell} {

            set CellAverageEdgeLength [$Dom getCellAverageEdgeLength $iCell]

            set minValue [expr {min($minValue, $CellAverageEdgeLength)}]
            set maxValue [expr {max($maxValue, $CellAverageEdgeLength)}]
          }
        }
      }

      set ToWrite "Max spacing on the TRex BC = "
      puts $ToWrite$maxValue
      set ToWrite "Min spacing on the TRex BC = "
      puts $ToWrite$minValue

      set ToWrite "Number of layers needed for min spacing to reach isotropy = "
      set minNOfLayers [expr {int(log($minValue/$BLSpacing)/log($GRBL))}]
      puts $ToWrite$minNOfLayers

      set ToWrite "Number of layers needed for max spacing to reach set wall distance = "
      set maxNLayers [expr {1- $BLHeight * (1-$GRBL) / $BLSpacing}]
      set maxNLayers [expr {int(log($maxNLayers)/log($GRBL))}]
      puts $ToWrite$maxNLayers

      set NMaxLayers [expr {max($minNOfLayers, $maxNLayers)}]

      $Volume setUnstructuredSolverAttribute TRexMaximumLayers $NMaxLayers
      $Volume setUnstructuredSolverAttribute TRexIsotropicHeight $IsotropicHeight
      $Volume setUnstructuredSolverAttribute TRexFullLayers $NFullLayers

    }

    if { $gridType == "Euler" } {
      $Volume setUnstructuredSolverAttribute TRexMaximumLayers 0
      $Volume setUnstructuredSolverAttribute TRexFullLayers 0
    }

    $Volume setSizeFieldDecay 0.8
    $Volume setUnstructuredSolverAttribute EdgeMaximumLength [lindex $MaxEdge $i]
    $_TMP(mode_1) end
    unset _TMP(mode_1)
    pw::Application markUndoLevel Solve
  }


  for {set i 0} {$i < [llength $Volumes]} {incr i} {
    set Volume [lindex $Volumes $i]
    set _TMP(mode_1) [pw::Application begin UnstructuredSolver [list $Volume]]
      $_TMP(mode_1) setStopWhenFullLayersNotMet true
      $_TMP(mode_1) setAllowIncomplete true

      set source [pw::SourceEntity getByName src-1]
      $Volume excludeSource $source false

      set ToWrite "Initializing block..."
      set ToWrite_1 [lindex $existingDomainsOrBlocks $i]
      set ToWrite_2 "..."
      puts $ToWrite$ToWrite_1$ToWrite_2

      $_TMP(mode_1) run Initialize

      set ToWrite "Smoothing block..."
      set ToWrite_1 [lindex $existingDomainsOrBlocks $i]
      set ToWrite_2 "..."
      puts $ToWrite$ToWrite_1$ToWrite_2

      $Volume setUnstructuredSolverAttribute WCNSmoothConvergenceCostThreshold 0.9
      $Volume setUnstructuredSolverAttribute WCNSmoothCostAngleThreshold 160
      $Volume setUnstructuredSolverAttribute WCNSmoothRelaxationFactor 0.15
      ###$_TMP(mode_1) run Smooth 400
    $_TMP(mode_1) end
    unset _TMP(mode_1)
    pw::Application markUndoLevel Initialize

    set ToWrite "Block created. Number of Points = "
    set NPoints [lindex [$Volume getDimensions] 0]
    set ToWrite_2 ". Number of Cells = "
    set NCells [$Volume getCellCount]
    puts $ToWrite$NPoints$ToWrite_2$NCells

  }
}


if { $OutputFormat == "SU2" } {

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

}

if { $OutputFormat == "CGNS" } {

  set Extension ".cgns"

  pw::Application setCAESolver CGNS $nDim
  pw::Application markUndoLevel {Select Solver}

  # pw::Application setCAESolverAttribute CGNS.FileType adf
  pw::Application setCAESolverAttribute CGNS.FileType hdf5
  # pw::Application setCAESolverAttribute CGNS.Units Meters
  # pw::Application setCAESolverAttribute CGNS.Units Centimeters
  # pw::Application setCAESolverAttribute CGNS.Units Millimeters
  # pw::Application setCAESolverAttribute CGNS.Units Feet
  # pw::Application setCAESolverAttribute CGNS.Units Inches
  pw::Application setCAESolverAttribute CGNS.Units UserDefined
  pw::Application setCAESolverAttribute CGNS.UseFamilyConditions false

  # set Extension ".pw"
  # pw::Application save $path2Mesh$RefLevel$Extension

  if { $nDim == 2 } {
    set _TMP(mode_1) [pw::Application begin CaeExport [pw::Entity sort $Domains]]
      $_TMP(mode_1) initialize -strict -type CAE $path2Mesh$Extension
      $_TMP(mode_1) setAttribute FilePrecision Double
      $_TMP(mode_1) setAttribute UnstructuredInterface NodeToNode
      $_TMP(mode_1) verify
      $_TMP(mode_1) write
    $_TMP(mode_1) end
    unset _TMP(mode_1)
  } else {
    set _TMP(mode_1) [pw::Application begin CaeExport [pw::Entity sort $Volumes]]
      $_TMP(mode_1) initialize -strict -type CAE $path2Mesh$Extension
      $_TMP(mode_1) setAttribute FilePrecision Double
      $_TMP(mode_1) setAttribute UnstructuredInterface NodeToNode
      $_TMP(mode_1) verify
      $_TMP(mode_1) write
    $_TMP(mode_1) end
    unset _TMP(mode_1)
  }

}

set Extension ".pw"
pw::Application save $path2Mesh$Extension
