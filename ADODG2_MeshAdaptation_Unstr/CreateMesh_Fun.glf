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
