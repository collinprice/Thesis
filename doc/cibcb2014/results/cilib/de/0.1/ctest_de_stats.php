<?php

$files = array_filter(glob('ctest*/data/exafs-de-*.txt'));

$data = array();

foreach($files as $file) {
	
	//echo "$file\n";
	$content = file_get_contents($file);
	$rows = explode("\n",$content);
	
	$counter = 0;
	foreach($rows as $row) {
		
		if (empty($row) || $row[0] == "#") continue;
		
		if (!array_key_exists($counter, $data)) $data[$counter] = 0;
		$cols = explode(" ", $row);
		$data[$counter++] += $cols[1];
		//echo "$cols[1]\n";
	}
}

$contents = "";
for($i = 0; $i < count($data); ++$i) {
	$data[$i] /= count($files);
	$contents .= "$i,$data[$i]\n";
}

file_put_contents("exafs-de-sum.csv",$contents);
//print_r($data);

