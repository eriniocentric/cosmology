use PDL;
use PGPLOT;
use PDL::NiceSlice;
use PDL::Graphics::PGPLOT::Window;
use PDL::Graphics::PGPLOTOptions ('set_pgplot_options');
set_pgplot_options('CharSize' => 1.4,'HardCH'=> 2.,'HardLW'=>3,'AspectRatio'=>1);

$OmegaM = 0.3; $H0 = 70;

$trials = 100;#number of bootstraps
($z_cat,$V_cat) = rcols("$ENV{HOME}/Software/Catalogs/comovingvolume.txt",0,1);
#$V_cat *= 90.75/121.0;
$V_cat *= 554.0/121.0;

($id,$spflag,$spectralclass,$z,$massKarl,$massKarl_err,$K,$Conf,$weight,$sfr2000,$sfrOII,$restUB,$gini,$assym,$fac,$facburst,$burstratio,$burstratio_avg,$burstratio_err,$bbtemp,$temp_avg,$temp_err,$massVIzK,$massVIzKerr,$zmaxVIzK,$massVIzK12,$massVIzK12err,$zmaxVIzK12,$massVIzK1234,$massVIzK1234err,$zmaxVIzK1234,$massVIzK1234_cce,$massVIzK1234_cceerr,$zmaxVIzK1234_cce,$rest3micronflux,$rest3micronflux_err,$massVIzK_nb,$massVIzK_nberr,$massVIzK12_nb,$massVIzK12_nberr,$massVIzK12_nb,$massVIzK12_nberr,$t,$t_err,$tburst,$tburst_err,$stellarburstratio,$stellarburstratio_err,$blendflag,$IRACflag) = rcols("$ENV{HOME}/Software/IRexcess/irexcess_gc.txt");

#combine SFRs, take average of two values if both exist

$sfr = -9.99 * ones(nelem($id));

$ix = which($sfr2000 == 0 & $sfrOII > 0);
$sfr($ix) .= $sfrOII($ix);

$ix = which($sfr2000 > 0 & $sfrOII == 0);
$sfr($ix) .= $sfr2000($ix);

$ix = which($sfr2000 > 0 & $sfrOII > 0);
$sfr($ix) .= ($sfr2000($ix) + $sfrOII($ix))/2;

$zlow  = pdl(0.8,1.0,1.3,1.6);
$zhigh = pdl(1.0,1.3,1.6,2.0);
$zplot = pdl(0.95,1.2,1.4,1.75);
$zerr  = pdl(0.15,0.1,0.1,0.25);

$nbin = nelem($zlow);

$massbins_s = pdl (9.0,10.2,10.8);
$massbins_e = pdl (10.2,10.8,11.5);
$nmsbin = nelem($massbins_s);

#$flag = which(($conf <= 1 | $zsp> 9) & $zsp>0.01);
#$z = pdl($zsp);
#$z->dice($flag) .= $zph->dice($flag);

$vol     = zeroes(nelem($id));
$sfrd     = zeroes($nmsbin,$nbin);
$sfrderr     = zeroes($nmsbin,$nbin);
$sfrd_bs     = zeroes($trials,$nbin);

$volnew     = zeroes(nelem($id));
$sfrdnew     = zeroes($nmsbin,$nbin);
$sfrdnewerr     = zeroes($nmsbin,$nbin);
$sfrdnew_bs     = zeroes($trials,$nbin);

$msbin = 0;

for ($msbin==0; $msbin < $nmsbin; $msbin++) {  
  
  for ($n=0;$n<=($trials-1);$n++) {
    $rnd = floor(random(nelem($id))*(nelem($id)));
    $rnd = sequence(nelem($id)) if ($n == 0);
    $wght        = $weight->dice($rnd);
    $ms      = $massVIzK->dice($rnd);
    $msnew  = $massVIzK1234_cce->dice($rnd);
    $zmx      = $zmaxVIzK->dice($rnd);
    $zmxnew = $zmaxVIzK1234_cce->dice($rnd);
    $sfr  = $sfr->dice($rnd);
    
    for ($bin = 0; $bin<=($nbin-1);$bin++) {
	for ($i=0;$i<=(nelem($id)-1);$i++) {
	    unless ($z(($i)) >= $zlow->(($bin)) & $z(($i)) <= $zhigh->(($bin))) {
	      $vol(($i)) .= 0;
	      $volnew(($i)) .= 0;
	      next;
	    }
	    
	    $z_l = $zlow(($bin));
	    $z_h = $zhigh(($bin));
	    
	    $z_h = $zmx(($i)) if ($zmx(($i)) > $z_l & $zmx(($i)) < $z_h);
	    $vol(($i)) .= $V_cat->(which(abs($z_h-$z_cat) < 0.005)->index(0))-$V_cat->index(which(abs($z_l-$z_cat) < 0.005)->index(0));
	    	    
	    $z_h = $zmxnew(($i)) if ($zmxnew(($i)) > $z_l & $zmxnew(($i)) < $z_h);
	    $volnew(($i)) .= $V_cat->(which(abs($z_h-$z_cat) < 0.005)->index(0))-$V_cat->index(which(abs($z_l-$z_cat) < 0.005)->index(0));
	    
	}#closes loop assigning volumes to each object
	
	$idx = which($vol > 1 & $ms > $massbins_s($msbin) & $ms < $massbins_e($msbin) & $K < 20.6 & $sfr > 0);       
	$sfrd_bs($n,$bin) .= sum(($sfr($idx))/(($wght($idx)*$vol($idx))));
	
	$idx = which($volnew > 1 & $msnew > $massbins_s($msbin) & $msnew < $massbins_e($msbin) & $K < 20.6 & $sfr >0);  
	$sfrdnew_bs($n,$bin) .=sum(($sfr($idx))/(($wght($idx)*$volnew($idx))));

      }#closes loop over bins
 
  }#closes loop over bootstrapping runs
  ($avg,$dum1,$dum2,$dum3,$dum4,$err) = statsover($sfrd_bs);
  $sfrd(($msbin),:) .= $avg;
  $sfrderr(($msbin),:) .= $err;

  ($avg,$dum1,$dum2,$dum3,$dum4,$err) = statsover($sfrdnew_bs);
  $sfrdnew(($msbin),:) .= $avg;
  $sfrdnewerr(($msbin),:) .= $err;

}#close loop over massbins

$logsfrderr = 0.4343 * $sfrderr / $sfrd;
$logsfrdnewerr = 0.4343 * $sfrdnewerr / $sfrdnew;

#dev '/xs';
dev "$ENV{HOME}/Software/Figures/sfrd.ps/vcps";
env 0.5,2,-3.5,-0.3,{Axis=>'LogY',YTitle=>'SFRD',XTitle=>'z'};

$msbin =0;
for ($msbin == 0; $msbin < $nmsbin; $msbin++) {  
  errb $zplot,log10($sfrd($msbin,:)),$zerr,$logsfrderr($msbin,:),{Symbol=>17};
}

#plot total SFRD

points $zplot,log10(sumover($sfrd)),{Symbol=>18,Color=>yellow};

close_window;
