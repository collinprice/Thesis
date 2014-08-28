<?php

$file = "results.csv";
$n = 180;

$content = file_get_contents($file);
$rows = explode("\n",$content);

$new_content = "";
$counter = 0;
foreach($rows as $row) {
        
        $cols = explode(",", $row);
        if (count($cols) != 3) continue;

        $cols[1] = sqrt($cols[1]/$n);
        $cols[2] = sqrt($cols[2]/$n);

        $new_content .= "$cols[0],$cols[1],$cols[2]\n";
}

file_put_contents("results-fixed-2.csv",$new_content);
