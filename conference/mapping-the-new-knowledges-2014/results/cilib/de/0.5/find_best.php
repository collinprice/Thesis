<?php

$files = array_filter(glob('test*/run*/results.csv'));

$best = 9999;

foreach($files as $file) {
	
	$line = `tail -n 1 $file`;

	$parts = explode(",", $line);
	$score = $parts[0];
	
	if ($score < $best) {
		$best = $score;
		$best_file = $file;
	}
}

echo $best_file . " " . $best . "\n";
