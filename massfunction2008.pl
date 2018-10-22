use PDL;
use PGPLOT;
use PDL::NiceSlice;
use PDL::Graphics::PGPLOT::Window;
use PDL::Graphics::PGPLOTOptions ('set_pgplot_options');
set_pgplot_options('CharSize' => 1.4,'HardCH'=> 2.,'HardLW'=>3,'AspectRatio'=>1);

$OmegaM = 0.3; $H0 = 70;

$trials = 10;#number of bootstraps
($z_cat,$V_cat) = rcols("$ENV{HOME}/Software/Catalogs/comovingvolume.txt",0,1);
#$V_cat *= 90.75/121.0;
$V_cat *= 554.0/121.0;

#($id,$spflag,$spectralclass,$z,$massKarl,$massKarl_err, $K, $Conf,$weight,$sfr2000,$sfrOII,$restUB,$gini,$assym,$fac,$facburst,$burstratio,$burstratio_avg,$burstratio_err,$bbtemp,$temp_avg,$temp_err,$massVIzK,$massVIzKerr,$zmaxVIzK,$massVIzK12,$massVIzK12err,$zmaxVIzK12,$massVIzK1234,$massVIzK1234err,$zmaxVIzK1234,$massVIzK1234_cce,$massVIzK1234_cceerr,$zmaxVIzK1234_cce,$rest4micronflux,$rest4micronflux_err) = rcols("$ENV{HOME}/Software/IRexcess/irexcess_cce.txt");

($id,$spflag,$spectralclass,$z,$massKarl,$massKarl_err,$K,$Conf,$weight,$sfr2000,$sfrOII,$restUB,$gini,$assym,$fac,$facburst,$burstratio,$burstratio_avg,$burstratio_err,$bbtemp,$temp_avg,$temp_err,$massVIzK,$massVIzKerr,$zmaxVIzK,$massVIzK12,$massVIzK12err,$zmaxVIzK12,$massVIzK1234,$massVIzK1234err,$zmaxVIzK1234,$massVIzK1234_cce,$massVIzK1234_cceerr,$zmaxVIzK1234_cce,$rest3micronflux,$rest3micronflux_err,$massVIzK_nb,$massVIzK_nberr,$massVIzK12_nb,$massVIzK12_nberr,$massVIzK12_nb,$massVIzK12_nberr,$t,$t_err,$tburst,$tburst_err,$stellarburstratio,$stellarburstratio_err,$blendflag,$IRACflag,$sfr2000sed,$sfrOIIsed,$sfr2000obs,$sfrOIIobs) = rcols("$ENV{HOME}/Software/IRexcess/irexcess_gc.txt");

($id,$name) = rcols("$ENV{HOME}/Software/IRexcess/irexcess_gc_names.txt");
#$weight       = (554.0/90.75)*$weight;#convert from LCIRS to GDDS
#$weight = (121.0/554.0)*$weight; 
#
#$zlow  = pdl(0.8,1.0,1.3,1.6);#erin's zbins
#$zhigh = pdl(1.0,1.3,1.6,2.0);
#$zplot = pdl(0.95,1.2,1.4,1.75);
#$zerr  = pdl(0.15,0.1,0.1,0.25);
#$completenesslimits = pdl(9.81,10.08,10.35,10.59);
#$who = "erin";

#$zlow  = pdl(0.8,1.1,1.4,1.7);#uses bob's zbins
#$zhigh = pdl(1.1,1.4,1.7,2.0);
#$zplot = pdl(0.95,1.25,1.55,1.85);
#$zerr  = pdl(0.15,0.15,0.15,0.15);
#$who = "bob";
#$completenesslimits = pdl(9.9,10.17,10.43,10.59);
$zpgplot = pdl(1.0,1.3,1.6,2.0);
$zfonplot = pdl(0.8,1.0,1.3,1.6);

#$zlow  = pdl(0.8,1.1,1.4);#uses bob's zbins
#$zhigh = pdl(1.1,1.4,1.7);
#$zplot = pdl(0.95,1.25,1.7);
#$zerr  = pdl(0.15,0.15,0.3);
#$who = "bob2";

#$zlow  = pdl(0.8,1.1,1.3,1.6);#uses paper 3 zbins
#$zhigh = pdl(1.1,1.3,1.6,2.0);
#$zplot = pdl(0.95,1.2,1.45,1.8);
#$zerr  = pdl(0.15,0.1,0.15,0.2);
#$who = "karl";

#$zlow  = pdl(0.8,1.2,1.6);
#$zhigh = pdl(1.2,1.6,2.0);
#$zplot = pdl(1.0,1.4,1.8);
#$zerr  = pdl(0.2,0.2,0.2);
#$who = "erin2";

$zlow  = pdl(0.7,1.0,1.5);
$zhigh = pdl(1.0,1.5,2.0);
$zplot = pdl(0.85,1.25,1.75);
$zerr  = pdl(0.15,0.25,0.25);
$who = "k20";
$completenesslimits = pdl(9.9,10.4,10.8);

$nbin = nelem($zlow);

#$massbins_low  = pdl( 8.4,8.8,9.2,9.6,10.0,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11.0,11.2);
#$massbins_high = pdl( 8.8,9.2,9.6,10.0,10.1,10.2,10.3,10.4,10.5,10.6,10.7,10.8,10.9,11.0,11.2,11.6);
#$massbins_low = pdl( 8.4,8.7,9.0,9.3,9.6,9.9, 10.2,10.4,10.6,10.8,11.0,11.2);
#$massbins_high = pdl(8.7,9.0,9.3,9.6,9.9,10.2,10.4,10.6,10.8,11.0,11.2,11.6);

$massbins_low = sequence(11)/4.+9.0;
$massbins_high = $massbins_low + 0.25;

$masserr = ($massbins_high-$massbins_low)/2.;
$massbins = ($massbins_low + $massbins_high)/2.;
$nmsbin = nelem($massbins);

#$flag = which(($conf <= 1 | $zsp> 9) & $zsp>0.01);
#$z = pdl($zsp);
#$z->dice($flag) .= $zph->dice($flag);

#$volKarl     = zeroes(nelem($idx));
$volVIzK     = zeroes(nelem($id));
$volVIzK12    = zeroes(nelem($id));
$volVIzK1234  = zeroes(nelem($id));
$volVIzK1234_cce = zeroes(nelem($id));

#$rhoKarl     = zeroes($nbin);
$rhoVIzK     = zeroes($nbin);
$rhoVIzK12     = zeroes($nbin);
$rhoVIzK1234  = zeroes($nbin);
$rhoVIzK1234_cce = zeroes($nbin);

#$rhoKarlavg     = zeroes($nbin);
$rhoVIzKavg     = zeroes($nbin);
$rhoVIzK12avg     = zeroes($nbin);
$rhoVIzK1234avg  = zeroes($nbin);
$rhoVIzK1234_cceavg = zeroes($nbin);

#$rhoKarlerr     = zeroes($nbin);
$rhoVIzKerr     = zeroes($nbin);
$rhoVIzK12err     = zeroes($nbin);
$rhoVIzK1234err  = zeroes($nbin);
$rhoVIzK1234_cceerr = zeroes($nbin);

#$rho_bsKarl     = zeroes($trials,$nbin);
$rho_bsVIzK     = zeroes($trials,$nbin);
$rho_bsVIzK12     = zeroes($trials,$nbin);
$rho_bsVIzK1234  = zeroes($trials,$nbin);
$rho_bsVIzK1234_cce = zeroes($trials,$nbin);

$phi_VIzK = zeroes($nmsbin,$nbin);
$phi_VIzK12 = zeroes($nmsbin,$nbin);
$phi_VIzK1234 = zeroes($nmsbin,$nbin);
$phi_VIzK1234_cce = zeroes($nmsbin,$nbin);

$phi_VIzK_err = zeroes($nmsbin,$nbin);
$phi_VIzK12_err = zeroes($nmsbin,$nbin);
$phi_VIzK1234_err = zeroes($nmsbin,$nbin);
$phi_VIzK1234_cce_err = zeroes($nmsbin,$nbin);

for ($n=0;$n<=($trials-1);$n++) {
    $rnd = floor(random(nelem($id))*(nelem($id)));
    $rnd = sequence(nelem($id)) if ($n == 0);
    $wght        = $weight->dice($rnd);
    $msVIzK      = $massVIzK->dice($rnd);
    $msVIzK12      = $massVIzK12->dice($rnd);
    $msVIzK1234   = $massVIzK1234->dice($rnd);
    $msVIzK1234_cce  = $massVIzK1234_cce->dice($rnd);
 #   $msKarl      = $massKarl->dice($rnd);
 #   $zmxKarl     = $zmaxKarl->dice($rnd);
    $zmxVIzK      = $zmaxVIzK->dice($rnd);
    $zmxVIzK12    = $zmaxVIzK12->dice($rnd);
    $zmxVIzK1234  = $zmaxVIzK1234->dice($rnd);
    $zmxVIzK1234_cce = $zmaxVIzK1234_cce->dice($rnd);

    for ($bin = 0; $bin<=($nbin-1);$bin++) {
	for ($i=0;$i<=(nelem($id)-1);$i++) {
	    unless ($z(($i)) >= $zlow->(($bin)) & $z(($i)) <= $zhigh->(($bin))) {
		#$volKarl(($i)) .= 0;
		$volVIzK(($i)) .= 0;
		$volVIzK12(($i)) .= 0;
		$volVIzK1234(($i)) .= 0;
		$volVIzK1234_cce(($i)) .= 0;
		next;
	    }
	    #$z_l = $zlow(($bin));
	    #$z_h = $zhigh(($bin));
	    #$z_h = $zmxKarl(($i)) if ($zmxKarl(($i)) > $z_l & $zmxKarl(($i)) < $z_h);
	    #$volKarl(($i)) .= $V_cat->(which(abs($z_h-$z_cat) < 0.005)->index(0))-$V_cat->index(which(abs($z_l-$z_cat) < 0.005)->index(0));	
	    
	    $z_l = $zlow(($bin));
	    $z_h = $zhigh(($bin));
	    
	    $z_h = $zmxVIzK(($i)) if ($zmxVIzK(($i)) > $z_l & $zmxVIzK(($i)) < $z_h);
	    $volVIzK(($i)) .= $V_cat->(which(abs($z_h-$z_cat) < 0.005)->index(0))-$V_cat->index(which(abs($z_l-$z_cat) < 0.005)->index(0));
	    
	    $z_h = $zmxVIzK12(($i)) if ($zmxVIzK12(($i)) > $z_l & $zmxVIzK12(($i)) < $z_h);
	    $volVIzK12(($i)) .= $V_cat->(which(abs($z_h-$z_cat) < 0.005)->index(0))-$V_cat->index(which(abs($z_l-$z_cat) < 0.005)->index(0));
	    
	    $z_h = $zmxVIzK1234(($i)) if ($zmxVIzK1234(($i)) > $z_l & $zmxVIzK1234(($i)) < $z_h);
	    $volVIzK1234(($i)) .= $V_cat->(which(abs($z_h-$z_cat) < 0.005)->index(0))-$V_cat->index(which(abs($z_l-$z_cat) < 0.005)->index(0));
	    
	    $z_h = $zmxVIzK1234_cce(($i)) if ($zmxVIzK1234_cce(($i)) > $z_l & $zmxVIzK1234_cce(($i)) < $z_h);
	    $volVIzK1234_cce(($i)) .= $V_cat->(which(abs($z_h-$z_cat) < 0.005)->index(0))-$V_cat->index(which(abs($z_l-$z_cat) < 0.005)->index(0));
	}#closes loop assigning volumes to each object
	
#	$id = which($volKarl > 1 & $msKarl > 10.5);
#	$rho_bsKarl($n,$bin) .= sum((10**$msKarl($id))/(($wght($id)*$volKarl($id))));
#	if ($n == 0) {
#	    $rhoKarl($bin) .= sum((10**$msKarl($id))/(($wght($id)*$volKarl($id))));
#	}
	
	$idx = which($volVIzK > 1 & $msVIzK > 10.5 & $K < 20.6);
       
	$rho_bsVIzK($n,$bin) .= sum((10**$msVIzK($idx))/(($wght($idx)*$volVIzK($idx))));
	
	if ($n == 0) {
	    $rhoVIzK($bin) .= sum((10**$msVIzK($idx))/(($wght($idx)*$volVIzK($idx))));
	}
	
	$idx = which($volVIzK12 > 1 & $msVIzK12 > 10.5 & $K < 20.6);
       
	$rho_bsVIzK12($n,$bin) .= sum((10**$msVIzK12($idx))/(($wght($idx)*$volVIzK12($idx))));
	
	if ($n == 0) {
	    $rhoVIzK12($bin) .= sum((10**$msVIzK12($idx))/(($wght($idx)*$volVIzK12($idx))));
	}

	$idx = which($volVIzK1234 > 1 & $msVIzK1234 > 10.5 & $K < 20.6);
	$rho_bsVIzK1234($n,$bin) .= sum((10**$msVIzK1234($idx))/(($wght($idx)*$volVIzK1234($idx))));
	
	if ($n == 0) {
	    $rhoVIzK1234($bin) .= sum((10**$msVIzK1234($idx))/(($wght($idx)*$volVIzK1234($idx))));
	}
	
	$idx = which($volVIzK1234_cce > 1 & $msVIzK1234_cce > 10.5 & $K < 20.6);
	$rho_bsVIzK1234_cce($n,$bin) .= sum((10**$msVIzK1234_cce($idx))/(($wght($idx)*$volVIzK1234_cce($idx))));
	
	if ($n == 0) {
	    $rhoVIzK1234_cce($bin) .= sum((10**$msVIzK1234_cce($idx))/(($wght($idx)*$volVIzK1234_cce($idx))));
	}	

	#calculate mass density function per stellar mass bin
	
	if ($n == 0) {

	    for ($msbin = 0; $msbin < $nmsbin; $msbin++) {
		
		$idx = which( $K <= 20.6 & $volVIzK > 1 & $msVIzK > $massbins_low($msbin) & $msVIzK < $massbins_high($msbin));
		
		$phi_VIzK($msbin,$bin) .= sum(1/(($wght($idx)*$volVIzK($idx))))/0.25 + 1e-9;
		$phi_VIzK_err($msbin,$bin) .= sqrt(sum(1/(($wght($idx)*$volVIzK($idx)))**2))/0.25;
	    	if ($massbins_low($msbin) == 10.25) {
		    wcols "%6s %8.3f %8.3f %8.3f %8.3f", $id($idx), $msVIzK($idx), $zmxVIzK($idx), $wght($idx), $volVIzK($idx), "$ENV{HOME}/Software/massfunction/mscrosscheck.dat";
		}
		$idx = which( $K <= 20.6 & $volVIzK12 > 1 & $msVIzK12 > $massbins_low($msbin) & $msVIzK12 < $massbins_high($msbin));
		
		$phi_VIzK12($msbin,$bin) .= sum(1/(($wght($idx)*$volVIzK($idx))))/0.25 + 1e-9;
		$phi_VIzK12_err($msbin,$bin) .= sqrt(sum(1/(($wght($idx)*$volVIzK($idx)))**2))/0.25;
	
		$idx = which( $K <= 20.6 & $volVIzK1234 > 1 & $msVIzK1234 > $massbins_low($msbin) & $msVIzK1234 < $massbins_high($msbin));
		
		$phi_VIzK1234($msbin,$bin) .= sum(1/(($wght($idx)*$volVIzK1234($idx))))/0.25 + 1e-9;
		$phi_VIzK1234_err($msbin,$bin) .= sqrt(sum(1/(($wght($idx)*$volVIzK1234($idx)))**2))/0.25;

		$idx = which( $K <= 20.6 & $volVIzK1234_cce > 1 & $msVIzK1234_cce > $massbins_low($msbin) & $msVIzK1234_cce < $massbins_high($msbin));

		print "For zbin = $bin and massbinhigh = $msbin, objectindex = $idx\n";
		$phi_VIzK1234_cce($msbin,$bin) .= sum(1/(($wght($idx)*$volVIzK1234_cce($idx))))/0.25 + 1e-9;
		$phi_VIzK1234_cce_err($msbin,$bin) .= sqrt(sum(1/(($wght($idx)*$volVIzK1234_cce($idx)))**2))/0.25;

	    }#end of loop over massbins

	}#end of if statement to calculate mass function on first monte carlo loop

    }#closes loop over bins

}#closes loop over bootstrapping runs

#($rhoKarlavg,$dum1,$dum3,$dum4,$dum4,$rhoKarlerr) = statsover($rho_bsKarl);
($rhoVIzKavg,$dum1,$dum3,$dum4,$dum4,$rhoVIzKerr) = statsover($rho_bsVIzK);
($rhoVIzK12avg,$dum1,$dum3,$dum4,$dum4,$rhoVIzK12err) = statsover($rho_bsVIzK12);
($rhoVIzK1234avg,$dum1,$dum3,$dum4,$dum4,$rhoVIzK1234err) = statsover($rho_bsVIzK1234);
($rhoVIzK1234_cceavg,$dum1,$dum3,$dum4,$dum4,$rhoVIzK1234_cceerr) = statsover($rho_bsVIzK1234_cce);

#$win = PDL::Graphics::PGPLOT::Window->new(Device=>'/xs');
$win = PDL::Graphics::PGPLOT::Window->new(Device=>"$ENV{HOME}/Software/Figures/new2010/massdensity".$who."bins.ps/vcps",AspectRatio=>1);
$win->env(0,5.1,6.9,8.5,{Charsize=>1.7,AxisColour=>'Black',XTitle=>'Redshift',YTitle=>'log\d10 \u\gr\d*\u (M>10\u10.5\dM\d\(2281)\u)(M\d\(2281)\u Mpc\u-3\d)'});
$ix = sequence(4);
$yhigherr = log10($rhoVIzK+$rhoVIzKerr)-log10($rhoVIzK);
$ylowerr  = log10($rhoVIzK)-log10($rhoVIzK-$rhoVIzKerr);
$ix = which($ylowerr > 0);
$win->points($zplot($ix),log10($rhoVIzK($ix)),{SymbolSize=>1.5,Symbol=>16,Colour=>Black,Plotline=>1});
$win->errb($zplot($ix),log10($rhoVIzK($ix)),$zerr($ix),$zerr($ix),$ylowerr($ix),$yhigherr($ix),{Colour=>Black});


$yhigherr = log10($rhoVIzK12+$rhoVIzK12err)-log10($rhoVIzK12);
$ylowerr  = log10($rhoVIzK12)-log10($rhoVIzK12-$rhoVIzK12err);
$ix = which($ylowerr > 0);
$win->points($zplot($ix),log10($rhoVIzK12($ix)),{SymbolSize=>1.5,Symbol=>13,Colour=>Green,Plotline=>3});
$win->errb($zplot($ix),log10($rhoVIzK12($ix)),$zerr($ix),$zerr($ix),$ylowerr($ix),$yhigherr($ix),{Colour=>Green});


$yhigherr = log10($rhoVIzK1234+$rhoVIzK1234_err)-log10($rhoVIzK1234);
$ylowerr  = log10($rhoVIzK1234)-log10($rhoVIzK1234-$rhoVIzK1234_err);
$yhigherr = log10($rhoVIzK1234_cce+$rhoVIzK1234_cceerr)-log10($rhoVIzK1234_cce);
$ylowerr  = log10($rhoVIzK1234_cce)-log10($rhoVIzK1234_cce-$rhoVIzK1234_cceerr);
$ix = which($ylowerr > 0);
$win->points($zplot($ix),log10($rhoVIzK1234($ix)),{SymbolSize=>1.5,Symbol=>18,Colour=>blue,Plotline=>1,LineStyle=>5});
$win->errb($zplot($ix),log10($rhoVIzK1234($ix)),$zerr($ix),$zerr($ix),$ylowerr($ix),$yhigherr($ix),{Colour=>blue});

$yhigherr = log10($rhoVIzK1234_cce+$rhoVIzK1234_cceerr)-log10($rhoVIzK1234_cce);
$ylowerr  = log10($rhoVIzK1234_cce)-log10($rhoVIzK1234_cce-$rhoVIzK1234_cceerr);
$ix = which($ylowerr > 0);
$win->points($zplot($ix),log10($rhoVIzK1234_cce($ix)),{SymbolSize=>1.5,Symbol=>17,Colour=>red,Plotline=>1,LineStyle=>2});
$win->errb($zplot($ix),log10($rhoVIzK1234_cce($ix)),$zerr($ix),$zerr($ix),$ylowerr($ix),$yhigherr($ix),{Colour=>red});

#plots massdensity points from elsner et al 2008, fontana et al 2006, 

($z_other,$logrho) = rcols( "$ENV{HOME}/Software/massfunction/otherguys/md_els08.txt");
$logrho += 2*log10(0.55);
$win->points($z_other,$logrho,{Col=>black,Symbol=>plus,Plotline=>1});
($z_other,$logrho) = rcols( "$ENV{HOME}/Software/massfunction/otherguys/md_f06.txt");
$logrho += 2*log10(0.55);
#$win->points($z_other,$logrho,{Col=>black,Symbol=>cross,Plotline=>1});
($z_other,$logrho) = rcols( "$ENV{HOME}/Software/massfunction/otherguys/md_per08.txt");
$logrho += 2*log10(0.55);
$win->points($z_other,$logrho,{Col=>black,Symbol=>3,Plotline=>1});

$win->legend(["VIz'K","VIz'K[3.6][4.5]","VIz'K[3.6][4.5][5.6][8.0]","VIz'K[3.6][4.5][5.6][8.0] (with NIR)"],0.5,8.4,{Colour=>[black,green,blue,red],Symbol=>[16,13,18,17],Charsize=>0.9,SymbolSize=>2.5,TextFraction=>0.7});
$win->legend(["VIz'K","VIz'K[3.6][4.5]","VIz'K[3.6][4.5][5.6][8.0]","VIz'K[3.6][4.5][5.6][8.0] (with NIR)"],0.5,8.4,{Colour=>[black,green,blue,red],LineStyle=>[1,3,5,2],Charsize=>0.9,TextFraction=>0.7});
$win->close();

#figure out mass completeness limits


for ($bin = 0; $bin < $nbin; $bin++) {
    
    $win = PDL::Graphics::PGPLOT::Window->new(Device=>"$ENV{HOME}/Software/Figures/new2010/massfunction".$who."_z=".$zlow(($bin))."-".$zhigh(($bin)).".ps/vcps");
    $win->env(8.2,12,-6.0,-1,{AxisColour=>'Black',XTitle=>'M\d*\u (M\d\(2281)\u)',Axis=>'LogXY',YTitle=>'\u\gf\d\u (Mpc\u-3\d)'});
   
    $logphi_err = 0.4343 * $phi_VIzK_err->dice_axis(1,(($bin))) / $phi_VIzK->dice_axis(1,(($bin)));
    $logphi = log10($phi_VIzK->dice_axis(1,(($bin))));
    $ix = which($logphi > -9); #only plot nonzero phi values so that plotline works well 
    $win->points($massbins($ix),$logphi($ix),{SymbolSize=>1.5,Symbol=>16,Colour=>Black,Plotline=>1});
    $win->errb($massbins($ix),$logphi($ix),$masserr($ix),$logphi_err($ix),{Colour=>Black});

    $logphi_err = 0.4343 * $phi_VIzK12_err->dice_axis(1,(($bin))) / $phi_VIzK12->dice_axis(1,(($bin)));
    $logphi = log10($phi_VIzK12->dice_axis(1,(($bin))));
    $ix = which($logphi > -9); #only plot nonzero phi values so that plotline works well 
    $win->points($massbins($ix),$logphi($ix),{SymbolSize=>1.5,Symbol=>13,Colour=>Green,Plotline=>3});
    $win->errb($massbins($ix),$logphi($ix),$masserr($ix),$logphi_err($ix),{Colour=>green});

    $logphi_err = 0.4343 * $phi_VIzK1234_err->dice_axis(1,(($bin))) / $phi_VIzK1234->dice_axis(1,(($bin)));
    $logphi = log10($phi_VIzK1234->dice_axis(1,(($bin))));
    $ix = which($logphi > -9); #only plot nonzero phi values so that plotline works well 
    $win->points($massbins($ix),$logphi($ix),{SymbolSize=>1.5,Symbol=>18,Colour=>blue,Plotline=>1,LineStyle=>5});
    $win->errb($massbins($ix),$logphi($ix),$masserr($ix),$logphi_err($ix),{Colour=>blue});

    $logphi_err = 0.4343 * $phi_VIzK1234_cce_err->dice_axis(1,(($bin))) / $phi_VIzK1234_cce->dice_axis(1,(($bin)));
    $logphi = log10($phi_VIzK1234_cce->dice_axis(1,(($bin))));
    $ix = which($logphi > -9); #only plot nonzero phi values so that plotline works well 
    $win->points($massbins($ix),$logphi($ix),{SymbolSize=>1.5,Symbol=>17,Colour=>red,Plotline=>1,LineStyle=>2});
    $win->errb($massbins($ix),$logphi($ix),$masserr($ix),$logphi_err($ix),{Colour=>red});
    $win->line(pdl($completenesslimits(($bin)),$completenesslimits(($bin))),pdl(-7,0),{Colour=>'black',LineStyle=>'Dashed'});

#plot points from Perez-Gonzalez et al 2008
($zPG,$logmassPG,$logphiPG,$logphierrlowPG,$logphierrhighPG) = rcols("$ENV{HOME}/Software/massfunction/otherguys/mf_per08.txt");

#convert from Salpeter 1955 IMF to B&G

$logmassPG += log10(0.55);
$logphiPG += log10(0.55);

$ipg = which($zPG == $zpgplot($bin));
$xerr = ones(nelem($ipg)) * 0.1;
$win->errb($logmassPG($ipg),$logphiPG($ipg),$xerr,$xerr,$logphierrlowPG($ipg),$logphierrhighPG($ipg),{SymbolSize=>1,Color=>orange,PlotLine=>1,Symbol=>17});

    
#plot points from Perez-Gonzalez et al 2008
($zfon,$logmassfon,$logphifon) = rcols("$ENV{HOME}/Software/massfunction/otherguys/mf_f06.txt");

#convert from Salpeter 1955 IMF to B&G

$logmassfon += log10(0.55);
$logphifon += log10(0.55);

$ifon = which($zfon == $zfonplot($bin));
$xerr = ones(nelem($ifon)) * 0.1;
$win->points($logmassfon($ifon),$logphifon($ifon),{Color=>purple,PlotLine=>1,SymbolSize=>1});

# Plot Bell et al Mass-fn (converted to H0=70)
#    $phistar = 0.0044; $Ms = 10.93; $alpha = -0.83;
# Correct to B&G IMF
#    $Ms += log10(0.79);
#    $m = sequence(20000)/1000;
##    $win->pgsci(4);
#    $win->line( $m, log10( $phistar *(lschechter($m, $Ms, $alpha))));
#    $win->pgltext('t',-2,0.2,'Bell et al. z=0',0.05);
#    $win->text('Bell et al. z=0',-2,0.2);
# Plot Cole et al Mass-fn (converted to H0=70) which is Salpeter
#    $phistar = 0.003087; $Ms = 11.16; $alpha = -1.18;
# Correct to B&G IMF
#    $Ms += log10(0.55);
#    $m = sequence(20000)/1000;
##    $win->pgsci(5);
#    $win->line( $m, log10( $phistar *(lschechter($m, $Ms, $alpha))));
#    $win->text('Cole et al. z=0',-4,0.2);

# Plot Cole et al Data points
    
#    ($m, $phi, $dphi) = rcols "$ENV{HOME}/Software/massfunction/Cole-Data-MF.dat", 0, 3,4;
#    $m -= 2*log10($H0/100);  $phi *= ($H0/100)**3;  $dphi *= ($H0/100)**3;
## Correct to B&G IMF
#    $m += log10(0.55);
#    $win->points( $m, log10($phi)); 
#    $win->errb( $m, log10($phi), zeroes($m), zeroes($m),log10($phi)-log10($phi-$dphi), log10($phi+$dphi)-log10($phi));

#DONT Plot Baugh models
    $count = 0;
    while (0) {
#    for $zb (0,1.1, 1.45, 1.8) {
	($mbin, $n) = rcols "$ENV{HOME}/Software/massfunction/BaughModels/New/ms_Gran2000_sal_z$zb.dat";
	$mbin -= 2*log10($H0/100);
	# Correct to B&G IMF
	$mbin += log10(0.55);
#	$win->pgsls(4); 
#	$win->pgsci(2+$count);
	$nn = $n * ($H0/100)**3;
	$ii = which $nn>0;
	$win->line( $mbin($ii), log10($nn($ii)));
#	$win->pgltext('t',-2-$count++,0.7,"Baugh model z=$zb",0.05);
	$win->text("Baugh model z=$zb",-2-$count++,0.7);
#	$win->pgsls(1);
    }
   
    $win->legend(["VIz'K","VIz'K[3.6][4.5]","VIz'K[3.6][4.5][5.6][8.0]","VIz'K[3.6][4.5][5.6][8.0] (with NIR)"],8.4,-5.1,{Colour=>[black,green,blue,red],Symbol=>[16,13,18,17],Charsize=>0.9,SymbolSize=>1.5,TextFraction=>0.7});
    $win->legend(["VIz'K","VIz'K[3.6][4.5]","VIz'K[3.6][4.5][5.6][8.0]","VIz'K[3.6][4.5][5.6][8.0] (with NIR)"],8.4,-5.1,{Colour=>[black,green,blue,red],LineStyle=>[1,3,5,2],Charsize=>0.9,TextFraction=>0.7});
    
    $win->close();	    
}


$win2 = PDL::Graphics::PGPLOT::Window->new(Device=>"$ENV{HOME}/Software/Figures/new2010/massfunctionall".$who."bins.ps/vcps");
$win2->env(8.2,12,-5.3,-2,{AxisColour=>'Black',XTitle=>'M\d*\u (M\d\(2281)\u)',Axis=>'LogXY',YTitle=>'\u\gf\d\u (Mpc\u-3\d)'});

$col = pdl('black','green','red','orange');

for ($bin = 0; $bin < $nbin; $bin++) {
    
    $logphi_err = 0.4343 * $phi_VIzK1234_cce_err->dice_axis(1,(($bin))) / $phi_VIzK1234_cce->dice_axis(1,(($bin)));
    $logphi = log10($phi_VIzK1234_cce->dice_axis(1,(($bin))));
    $ix = which($logphi > -9 & $massbins >= $completenesslimits(($bin))); #only plot nonzero phi values so that plotline works
    # well and choose those with points above mass completness level of 80%
    $col = 'black';
    $col = 'green' if ( $bin == 1);
    $col = 'red' if ( $bin == 2);

    $win2->points($massbins($ix),$logphi($ix),{SymbolSize=>1.5,Symbol=>17,Colour=>$col,Plotline=>1,LineStyle=>1});
    $win2->errb($massbins($ix),$logphi($ix),$masserr($ix),$logphi_err($ix),{Colour=>$col});

    #$win2->line(pdl($completenesslimits(($bin)),$completenesslimits(($bin))),pdl(-4,0),{Colour=>'black',LineStyle=>'Dashed'});
    $win2->legend(["z = 0.7-1","z = 1-1.5","z = 1.5-2"],8.5,-4.5,{Colour=>[black,green,red],Symbol=>[17,17,17],Charsize=>1.1,TextFraction=>0.7});
}

$win2->close();

#write phi's and rhos to file

$phi_VIzK_out = $phi_VIzK->clump(-1);
$phi_VIzK_err_out = $phi_VIzK_err->clump(-1);
$phi_VIzK12_out = $phi_VIzK12->clump(-1);
$phi_VIzK12_err_out = $phi_VIzK12_err->clump(-1);
$phi_VIzK1234_out = $phi_VIzK1234->clump(-1);
$phi_VIzK1234_err_out = $phi_VIzK1234_err->clump(-1);
$phi_VIzK1234_cce_out = $phi_VIzK1234_cce->clump(-1);
$phi_VIzK1234_cce_err_out = $phi_VIzK1234_cce_err->clump(-1);

$masslowout = pdl($massbins_low,$massbins_low,$massbins_low)->clump(-1);
$masshighout = pdl($massbins_high,$massbins_high,$massbins_high)->clump(-1);

$zlowout = pdl($zlow(0)*ones(nelem($massbins)),$zlow(1)*ones(nelem($massbins)),$zlow(2)*ones(nelem($massbins)))->clump(-1);
$zhighout = pdl($zhigh(0)*ones(nelem($massbins)),$zhigh(1)*ones(nelem($massbins)),$zhigh(2)*ones(nelem($massbins)))->clump(-1);

#write massfunction to file phis.dat
wcols '  %4.2f-%4.2f & %4.2f-%4.2f & %6.3e  $\pm$  %6.3e  & %6.3e  $\pm$ %6.3e   & %6.3e  $\pm$ %6.3e   & %6.3e  $\pm$ %6.3e \\\ ', $zlowout,$zhighout,$masslowout,$masshighout,$phi_VIzK_out,$phi_VIzK_err_out,$phi_VIzK12_out,$phi_VIzK12_err_out,$phi_VIzK1234_out,$phi_VIzK1234_err_out,$phi_VIzK1234_cce_out,$phi_VIzK1234_cce_err_out, "$ENV{HOME}/Documents/Thesis/chapters/ch2/tables/massfunction.tab";
#write rho density to file rho.dat
wcols '  %4.2f-%4.2f & %6.3e  $\pm$  %6.3e  & %6.3e  $\pm$ %6.3e   & %6.3e  $\pm$ %6.3e   & %6.3e  $\pm$ %6.3e \\\ ', $zlow, $zhigh, $rhoVIzK,$rhoVIzKerr,$rhoVIzK12,$rhoVIzK12err,$rhoVIzK1234,$rhoVIzK1234err,$rhoVIzK1234_cce,$rhoVIzK1234_cceerr, "$ENV{HOME}/Documents/Thesis/chapters/ch2/tables/massdensity.tab";

# Schecter function in number per log mass
sub lschechter{
  my($M,$Ms,$alpha) = @_;
  my $x = 10**($M-$Ms);
  return 2.302585 * $x**($alpha+1) * exp(-$x);
}
