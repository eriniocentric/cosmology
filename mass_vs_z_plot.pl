use PDL; 
use PDL::NiceSlice; 
use PDL::Graphics::PGPLOT::Window; 
use PDL::Graphics::PGPLOTOptions ('set_pgplot_options');
set_pgplot_options('CharSize' => 0.8,'HardCH'=> 1.0,'HardLW'=>2.5);
chdir "$ENV{HOME}/Software/Masses/";

system "sed s/SA// <GDDSVIzK1234-masses.dat> tmp.txt";
system "sed s/-// <tmp.txt> temp.dat";
($idVIzK1234,$zVIzK1234,$KVIzK1234,$massVIzK1234,$massVIzK1234_err) = rcols("$ENV{HOME}/Software/Masses/temp.dat",0,1,2,3,4);

system "sed s/SA// <GDDSVIzK1234/zmax.dat> tmp.txt";
system "sed s/-// <tmp.txt> temp.txt";
($idVIzK1234z,$zVIzK1234max) = rcols("temp.txt",0,1);

system "sed s/SA// <GDDSVIzK-masses.dat> tmp.txt";
system "sed s/-// <tmp.txt> temp.dat";
($idVIzK,$zVIzK,$KVIzK,$massVIzK,$massVIzK_err) = rcols("$ENV{HOME}/Software/Masses/temp.dat",0,1,2,3,4);

system "sed s/SA// <GDDSVIzK/zmax.dat> tmp.txt";
system "sed s/-// <tmp.txt> temp.txt";
($idVIzKz,$zVIzKmax) = rcols("temp.txt",0,1);

system "sed s/SA// <GDDSVIzK1234_starburst-masses_9kpc.dat> tmp.txt";
system "sed s/-// <tmp.txt> temp.dat";
($idstarburst,$zstarburst,$Kstarburst,$massstarburst,$massstarburst_err) = rcols("$ENV{HOME}/Software/Masses/temp.dat",0,1,2,3,4);

system "sed s/SA// <GDDSVIzK1234_starburst_9kpc/zmax.dat> tmp.txt";
system "sed s/-// <tmp.txt> temp.txt";
($idstarburstz,$zstarburstmax) = rcols("temp.txt",0,1);

system "sed s/SA// <~/Software/Catalogs/Rest-colors.dat> tmp.txt";
system "sed s/-// <tmp.txt> temp.txt";
($id,$UBrest) = rcols("temp.txt",0,2);

$ix = zeroes(nelem($idVIzK));
for ($i=0;$i<=nelem($idVIzK)-1;$i++) {
    if (any $id == $idVIzK(($i))) {
	$tmp = which($id == $idVIzK(($i))); 
	$ix(($i)) .= $tmp(0);
    }
    else {
	$ix(($i)) .= -1;
    }
}

$UB = zeroes(nelem($idVIzK));
for ($i=0;$i<=(nelem($idVIzK)-1);$i++) {
    if ($ix(($i)) == -1) {
	$UB(($i)) .= -9.99;
	next;
    }
    $UB(($i)) .= $UBrest->index($ix(($i)));
}

($zml,$m1,$m2) = rcols("$ENV{HOME}/Software/Masses/mass-limits.dat",0,1,2);

$K = $KVIzK;
$z = $zVIzK;
				
$idx1 = which($K < 19.0 );
$idx2 = which($K >= 19.0 );
$idx3 = which($K >= 19.8 );

$opt3 = {Device => 'VIzK1234_vs_z.ps/cps'};
$win3 = PDL::Graphics::PGPLOT::Window->new($opt3);
$win3->errb($z($idx1),$massVIzK1234($idx1),$massVIzK1234_err($idx1), {SymbolSize=>1.4,Colour=>red,Symbol=>17, XRange=> [0.3,2.1],YRange=>[8.3,12]});
$win3->hold;
$win3->errb($z($idx2),$massVIzK1234($idx2),$massVIzK1234_err($idx2), {SymbolSize=>1.4,Colour=>green,Symbol=>16});
$win3->errb($z($idx3),$massVIzK1234($idx3),$massVIzK1234_err($idx3), {SymbolSize=>1.4,Colour=>blue,Symbol=>18});
$win3->line($zml,$m1,{Colour=>Black,LineStyle=>1});
$win3->line($zml,$m2,{Colour=>Black,LineStyle=>2});
$win3->legend(['K<19','19<K<19.8','19.8<K<20.6'], 1.6,8.7, {Colour=>[red,green,blue],Symbol=>[17,16,18]});
$win3->label_axes('Redshift','VIzK + IRAC 3.6, 4.5, 5.8, 8.0 log\d10\u M  (M\d\(2281)\u)');
$win3->release; 
$win3->close();

$idxUB1 = which($UB > 0.1 );
$idxUB2 = which($UB <= 0.1 & $UB >= -0.1 );
$idxUB3 = which($UB < -0.1) ;

$opt6 = {Device => "$ENV{HOME}/Software/paper/figures/VIzK1234_vs_VIzK_UBcuts.ps/cps"};
$win6 = PDL::Graphics::PGPLOT::Window->new($opt6);
$win6->errb($massVIzK($idxUB1),$massVIzK1234($idxUB1),$massVIzK_err($idxUB1),$massVIzK1234_err($idxUB1), {SymbolSize=>1.5,Colour=>red, Symbol=>17, XRange=> [8.5,11.7],YRange=>[8.5,11.7]});
$win6->hold;
$win6->errb($massVIzK($idxUB2),$massVIzK1234($idxUB2),$massVIzK_err($idxUB2),$massVIzK1234_err($idxUB2), {SymbolSize=>1.5,Colour=>orange,Symbol=>16});
$win6->errb($massVIzK($idxUB3),$massVIzK1234($idxUB3),$massVIzK_err($idxUB3),$massVIzK1234_err($idxUB3), {SymbolSize=>1.5,Colour=>blue,Symbol=>18});
$win6->line(sequence(20),sequence(20),{colour=>red});
$win6->legend(['  0.1 < (U-B)\drest\u','-0.1 < (U-B)\drest\u < 0.1','-0.1 > (U-B)\drest\u'], 10.2,9.1,{Colour=>[red,orange,blue],Symbol=>[17,16,18],Symbolsize=>1.7});
$win6->label_axes('VIzK  log\d10\u M  (M\d\(2281)\u) ','VIzK + IRAC log\d10\u M  (M\d\(2281)\u)');
$win6->release; 
$win6->close();

$opt6 = {Device => "$ENV{HOME}/Software/paper/figures/VIzK1234_vs_VIzK_spectraltypes.ps/cps"};
$win6 = PDL::Graphics::PGPLOT::Window->new($opt6);
$win6->errb($massVIzK($idxUB1),$massVIzK1234($idxUB1),$massVIzK_err($idxUB1),$massVIzK1234_err($idxUB1), {SymbolSize=>1.5,Colour=>red, Symbol=>17, XRange=> [8.5,11.7],YRange=>[8.5,11.7]});
$win6->hold;
$win6->errb($massVIzK($idxUB2),$massVIzK1234($idxUB2),$massVIzK_err($idxUB2),$massVIzK1234_err($idxUB2), {SymbolSize=>1.5,Colour=>orange,Symbol=>16});
$win6->errb($massVIzK($idxUB3),$massVIzK1234($idxUB3),$massVIzK_err($idxUB3),$massVIzK1234_err($idxUB3), {SymbolSize=>1.5,Colour=>blue,Symbol=>18});
$win6->line(sequence(20),sequence(20),{colour=>red});
$win6->legend(['  0.1 < (U-B)\drest\u','-0.1 < (U-B)\drest\u < 0.1','-0.1 > (U-B)\drest\u'], 10.2,9.1,{Colour=>[red,orange,blue],Symbol=>[17,16,18],Symbolsize=>1.7});
$win6->label_axes('VIzK  log\d10\u M  (M\d\(2281)\u) ','VIzK + IRAC log\d10\u M  (M\d\(2281)\u)');
$win6->release; 
$win6->close();

$opt6 = {Device => 'VIzK1234starburst_vs_VIzK1234_UBcuts.ps/cps'};
$win6 = PDL::Graphics::PGPLOT::Window->new($opt6);
$win6->errb($massVIzK1234($idxUB1),$massstarburst($idxUB1),$massVIzK1234_err($idxUB1),$massstarburst_err($idxUB1), {SymbolSize=>1.4,Colour=>red, Symbol=>17, XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win6->hold;
$win6->errb($massVIzK1234($idxUB2),$massstarburst($idxUB2),$massVIzK1234_err($idxUB2),$massstarburst_err($idxUB2), {SymbolSize=>1.4,Colour=>green,Symbol=>16});
$win6->errb($massVIzK1234($idxUB3),$massstarburst($idxUB3),$massVIzK1234_err($idxUB3),$massstarburst_err($idxUB3), {SymbolSize=>1.4,Colour=>blue,Symbol=>18});
$win6->line(sequence(20),sequence(20),{colour=>red});
$win6->legend(['  0.1 < UB\drest\u','-0.1 < UB\drest\u < 0.1','-0.1 > UB\drest\u'], 10.5,9.3,{Colour=>[red,green,blue],Symbol=>[17,16,18]});
$win6->label_axes('VIzK + IRAC log\d10\u M  (M\d\(2281)\u)','VIzK + IRAC + StarBurst log\d10\u M  (M\d\(2281)\u)');
$win6->release; 
$win6->close();

$idz1 = which($z < 1.2 );
$idz2 = which($z < 1.4 & $z >= 1.2 );
$idz3 = which($z < 2.0 & $z >= 1.4 );

$opt7 = {Device => 'VIzK1234_vs_VIzK_zcuts.ps/cps'};
$win7 = PDL::Graphics::PGPLOT::Window->new($opt7);
$win7->points($massVIzK($idz1),$massVIzK1234($idz1), {SymbolSize=>1.4,Colour=>blue, Symbol=>18,XRange=> [9.0,11.7],YRange=>[9.0,11.7]});
$win7->hold;
$win7->points($massVIzK($idz2),$massVIzK1234($idz2), {SymbolSize=>1.4,Colour=>green,Symbol=>16});
$win7->points($massVIzK($idz3),$massVIzK1234($idz3), {SymbolSize=>1.4,Colour=>red,Symbol=>17});
$win7->line(sequence(20),sequence(20),{colour=>red});
$win7->legend(['0.5 > z > 1.2','1.2 > z > 1.4','1.4 > z > 2'], 10.5,9.3,{Colour=>[blue,green,red],Symbol=>[18,16,17]});
$win7->label_axes('VIzK Mass','VIzK1234 Log Mass');
$win7->release; 
$win7->close();

$massdif_UB1 =  $massVIzK1234($idxUB1) - $massVIzK($idxUB1);
$massdif_UB2 =  $massVIzK1234($idxUB2) - $massVIzK($idxUB2);
$massdif_UB3 =  $massVIzK1234($idxUB3) - $massVIzK($idxUB3);

($mean1,$d2,$d3,$d4,$d5,$d6,$rms1) = stats($massdif_UB1);
($mean2,$d2,$d3,$d4,$d5,$d6,$rms2) = stats($massdif_UB2);
($mean3,$d2,$d3,$d4,$d5,$d6,$rms3) = stats($massdif_UB3);

$massdif_all = $massVIzK1234-$massVIzK;

($meanall,$d2,$d3,$d4,$d5,$d6,$rmsall) = stats($massdif_all);

system "mv *.ps $ENV{HOME}/Software/Figures/xps";
