<?php

$file = "exafs-pso-sum.csv";
$n = 181;

$content = file_get_contents($file);
$rows = explode("\n",$content);

$new_content = "";
$counter = 0;
foreach($rows as $row) {
        
        $cols = explode(",", $row);
        if (count($cols) != 2) continue;

        $cols[1] = sqrt($cols[1]/$n);

        $new_content .= "$cols[0],$cols[1]\n";
}

file_put_contents("results-fixed.csv",$new_content);