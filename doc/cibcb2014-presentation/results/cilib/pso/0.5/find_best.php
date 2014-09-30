<?php

$files = array_filter(glob('ctest*/data/exafs-pso-*.txt'));

$best = 9999;

foreach($files as $file) {
	
	$line = `tail -n 1 $file`;

	$parts = explode(" ", $line);
	$score = $parts[1];
	
	if ($score < $best) {
		$best = $score;
		$best_file = $file;
	}
}

echo $best_file . " " . $best . "\n";
