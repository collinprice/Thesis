<?php

$population = [
	[2,4],
	[2,10],
	[3,4],
	[4,3],
	[5,10],
];

function arrayContainsArray(array $container, array $element) {

	for ($i=0; $i < count($container); $i++) { 
		
		if (arraysAreEqual($container[$i], $element)) return true;
	}

	return false;
}

function arraysAreEqual(array $a, array $b) {

	for ($i=0; $i < count($a); $i++) { 
		if ($a[$i] != $b[$i]) return false;
	}

	return true;
}

/*
	return true = A dominates B
	return false = B dominates A
	return null = same rank

*/
function A_dominates_B(array $A, array $B) {

	$at_least_one_better = false;
	$failed = false;

	for ($i=0; $i < count($A); $i++) { 
		
		if ($A[$i] < $B[$i]) {
			$at_least_one_better = true;
			continue;
		} else {
			$failed = true;
		}
	}

	if (!$failed) {
		return true;
	} elseif ($failed && $at_least_one_better) {
		return null;
	} else {
		return false;
	}
}

$ranks = array();

while(count($population) > 0) {

	$current_rank = array();

	for ($i=0; $i < count($population); $i++) { 
		
		$keep = true;

		for ($j=0; $j < count($population); $j++) { 
		
			if ($i == $j) continue;
			// print_r($population[$i]);
			// print_r($population[$j]);
			$result = A_dominates_B($population[$i], $population[$j]);

			if ($result === true) {
				// echo "Dominates\n";
			} elseif ($result === null) {
				// echo "Equal\n";
			} else {
				$keep = false;
				// echo "Dominated By\n";
			}
		}

		// echo $keep ? "KEEP\n" : "NOPE\n";
		// echo "\n";

		if ($keep) {
			$current_rank[] = $population[$i];
		}
		
	}

	$new_population = array();

	for ($i=0; $i < count($population); $i++) { 
		
		// Check if current rank contains this index
		if (!arrayContainsArray($current_rank, $population[$i])) {
			$new_population[] = $population[$i];
		}
	}
	$ranks[] = $current_rank;
	$population = $new_population;

	print_r($new_population);

	// break;
}

print_r($ranks);