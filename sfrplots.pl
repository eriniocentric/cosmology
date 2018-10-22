
use PDL;
use lib "$ENV{HOME}/Software";

use PDL::NiceSlice;
use PDL::Graphics::PGPLOT;
use PGPLOT;
use PDL::Graphics::PGPLOTOptions ('set_pgplot_options');
set_pgplot_options('CharSize' => 2.3,'HardCH'=> 1.5,'HardLW'=>3,'AspectRatio'=>0.7,'TightLabels'=>1);
use KGB::Cosmology;
use KGB::PegUtils;
use KGB::SpecUtils;
use KGB::Dust;
use PDL::Fit::Polynomial;

$Law = "SMC"; # Dust Law MW/SMC

($ID,$spflag,$SpType,$z,$Conf,$sfr2000,$sfrOII,$sfr2000mistake,$sfrOIImistake,$logM,$logMerr,$zmax,$rest3micronflux,$rest3micronflux_err,$IRACflag,$blendflag,$SpWeight,$AV,$AVerr) = rcols("$ENV{HOME}/Software/Catalogs/gdds-sfr.dat");

($z2,$logM2,$logM2err,$FOII,$FOIIerr,$FOIIflag,$sfrOII2,$sfr20002) = rcols("$ENV{HOME}/Software/Catalogs/gdds-sfr.dat",19,21,22,23,24,25,26,27); 
#properly correct SFRs with AV extinction from SED fits

$pc = 3.086E16; # One parsec in meters

$DL = lumdist($z);
       
$FOIIobs = $FOII * 4*$Pi * ($DL*1E6*$pc)**2;

$sfrHalpha = 

$attn = peidust($Law, 6565)/peidust($Law, 5500);
$AVHalpha = 1;

$sfrOIIobs = $sfrOII * ( 10**(-0.4 * $AVHalpha ) ); #convert to observed SFR (not extinction corrected)

$sfrOIIsed = $sfrOIIobs / (10**(-0.4 * $attn * 2. * $AV )); #Apply Av measured from SED fits to Halpha SFR

$AV2000 = 2.2;

$attn = peidust($Law, 2000)/peidust($Law, 5500);

$sfr2000obs = $sfr2000 * (10**(-0.4 * $AV2000));#convert to observed SFR (not extinction corrected)
$sfr2000sed = $sfr2000obs / (10**(-0.4 * $attn * $AV)); #Apply Av measured from SED fits to UV2000 SFR

$attn = peidust($Law,2000) / peidust($Law,1500); # use A_1500-logM relationship from Panella et al

$A_1500 = 4.07 * ($logM - log10(0.55) ) - 39.32; #convert logM from BG03 to Salpeter 
$sfr2000pan = $sfr2000obs / (10**(-0.4 * $attn * $A_1500));

$attn = peidust($Law,6565) / peidust($Law,1500); # use A_1500-logM relationship from Panella et al
$sfrOIIpan = $sfr2000obs / (10**(-0.4 * $attn * $A_1500));


$ix = which( $rest3micronflux > $rest3micronflux_err & $sfrOIIsed > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag == 0 );

$sfr = $sfrOIIsed($ix);
$r3mu = $rest3micronflux($ix);
$r3mu_err = $rest3micronflux_err($ix);
$r3mu_wt = 1 / (0.4343 * $r3mu_err / $r3mu);

#fit again for the 3 micron excess

$corrcoef_r3mu = sumover( ($sfr - average($sfr)) * ($r3mu-average($r3mu))) / sqrt( sumover(($sfr - average($sfr))**2) * sumover(($r3mu - average($r3mu))**2));

#monte carlo $br and $ssfr to get errors on fit

$mctrials = 100;
$r3muint_trials = zeroes($mctrials);
$r3muslope_trials = zeroes($mctrials);

for ($mc = 0; $mc < $mctrials; $mc++) {

    $rnd = floor(random(nelem($ix))*(nelem($ix)));
    $rnd = sequence(nelem($ix)) if ($mc == 0);
    $r3mu_mc       = $r3mu->dice($rnd);
    $wt_mc        = $r3mu_wt->dice($rnd);
    $sfr_mc       = $sfr->dice($rnd);

    ($dummy,$coefs_mc) = fitpoly1d(log10($sfr_mc),log10($r3mu_mc),2,{Weights=> $wt_mc});
    $r3muint_trials($mc) .= $coefs_mc(0);
    $r3muslope_trials($mc) .= $coefs_mc(1);
  }

$slope = average($r3muslope_trials);
$int = 10**(average($r3muint_trials));

$opt = {Device=> "$ENV{HOME}/Software/paper/figures/LNIR_vs_sfr.ps/vcps"};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->env(-1.3,2.3,29.5,32.8,{Axis=>'LogXY',XTitle=>'SFR (M\d\(2281)\u/yr)',YTitle=>'L\d3\gmm\u (Greybody/PAH) (W/\A)',AspectRatio=>1});

$Lsun = 3.826e26; #in Watts
#write flux in solar units
#$rest3micronflux /= $Lsun; #in Lsun
#$rest3micronflux_err /= $Lsun;#in Lsun
$dlogrest3micronflux_err = 0.4343 * $rest3micronflux_err/ $rest3micronflux;

#sptype cut

$idx1 = which( $rest3micronflux > $rest3micronflux_err & $spflag == 0 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);
$idx2 = which( $rest3micronflux > $rest3micronflux_err & $spflag == 1 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);
$idx3 = which( $rest3micronflux > $rest3micronflux_err & $spflag == 2 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);
$idx4 = which( $rest3micronflux > $rest3micronflux_err & $spflag == -1 & $sfr2000 > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);

$win->errb(log10($sfr2000($idx1)), log10($rest3micronflux($idx1)), undef, $dlogrest3micronflux_err($idx1), {Symbol=>17,Colour=>'blue',SymbolSize=>2});
$win->errb(log10($sfr2000($idx2)), log10($rest3micronflux($idx2)), undef, $dlogrest3micronflux_err($idx2), {Symbol=>17,Colour=>'green',SymbolSize=>2});
$win->errb(log10($sfr2000($idx3)), log10($rest3micronflux($idx3)), undef, $dlogrest3micronflux_err($idx3), {Symbol=>17,Colour=>'red',SymbolSize=>2});
$win->errb(log10($sfr2000($idx4)), log10($rest3micronflux($idx4)), undef, $dlogrest3micronflux_err($idx4), {Symbol=>17,Colour=>'orange',SymbolSize=>2});

$idx1 = which( $rest3micronflux > $rest3micronflux_err & $spflag == 0 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);
$idx2 = which( $rest3micronflux > $rest3micronflux_err & $spflag == 1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);
$idx3 = which( $rest3micronflux > $rest3micronflux_err & $spflag == 2 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);
$idx4 = which( $rest3micronflux > $rest3micronflux_err & $spflag == -1 & $sfrOII > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);

$win->errb(log10($sfrOII($idx1)), log10($rest3micronflux($idx1)), undef, $dlogrest3micronflux_err($idx1), {Symbol=>18,Colour=>'blue',SymbolSize=>2});
$win->errb(log10($sfrOII($idx2)), log10($rest3micronflux($idx2)), undef, $dlogrest3micronflux_err($idx2), {Symbol=>18,Colour=>'green',SymbolSize=>2});
$win->errb(log10($sfrOII($idx3)), log10($rest3micronflux($idx3)), undef, $dlogrest3micronflux_err($idx3), {Symbol=>18,Colour=>'red',SymbolSize=>2});
$win->errb(log10($sfrOII($idx4)), log10($rest3micronflux($idx4)), undef, $dlogrest3micronflux_err($idx4), {Symbol=>18,Colour=>'orange',SymbolSize=>2});

$UL3micronflux = log10($rest3micronflux_err) + $dlogrest3micronflux_err;

#plot fit

$sfr = 10**(pdl(-4,-3,-2,1,0,1,2,3)); 
$r3mu = $int * $sfr**($slope);

$win->line(log10($sfr),log10($r3mu));

#use order of magnitude approximation Lexcess in Lband = 0.642 * (0.15) * sfr * lifetime
$Lexcess = 350; #calculated in circumstellar/disksize.pl
$flux3micron = $Lsun / 3e4 * $sfr * 1e6 * $Lexcess; #in units of W/A
$win->line(log10($sfr),log10($flux3micron),{Color=>green,Linestyle=>solid});
$win->text("Toy CS Model",-0.9,29.6,{Angle=>30,Color=>green,Charsize=>1});
$win->legend(["Red Evolved","Intermediate","Star Forming","Mixed"],-1.0,32.5,{Colour=>['Red','Green','Blue','orange'],Symbol=>[17,17,17,17],Charsize=>1.1,SymbolSize=>1.5,TextShift=>-1.65});

$win->close();

$opt = {Device=> "$ENV{HOME}/Software/paper/figures/LNIR_vs_sfrex.ps/vcps"};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->env(-1.3,2.3,29.5,32.8,{Axis=>'LogXY',XTitle=>'SFR (M\d\(2281)\u/yr - extinction corrected)',YTitle=>'L\d3\gmm\u (Greybody/PAH) (W/\A)',AspectRatio=>1});

$Lsun = 3.826e26; #in Watts
#write flux in solar units
#$rest3micronflux /= $Lsun; #in Lsun
#$rest3micronflux_err /= $Lsun;#in Lsun
$dlogrest3micronflux_err = 0.4343 * $rest3micronflux_err/ $rest3micronflux;

#sptype cut


$idx1 = which( $rest3micronflux > $rest3micronflux_err & $spflag == 0 & $sfr2000pan > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);
$idx2 = which( $rest3micronflux > $rest3micronflux_err & $spflag == 1 & $sfr2000pan > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);
$idx3 = which( $rest3micronflux > $rest3micronflux_err & $spflag == 2 & $sfr2000pan > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);
$idx4 = which( $rest3micronflux > $rest3micronflux_err & $spflag == -1 & $sfr2000pan > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);

$win->errb(log10($sfr2000pan($idx1)), log10($rest3micronflux($idx1)), undef, $dlogrest3micronflux_err($idx1), {Symbol=>17,Colour=>'blue',SymbolSize=>2});
$win->errb(log10($sfr2000pan($idx2)), log10($rest3micronflux($idx2)), undef, $dlogrest3micronflux_err($idx2), {Symbol=>17,Colour=>'green',SymbolSize=>2});
$win->errb(log10($sfr2000pan($idx3)), log10($rest3micronflux($idx3)), undef, $dlogrest3micronflux_err($idx3), {Symbol=>17,Colour=>'red',SymbolSize=>2});
$win->errb(log10($sfr2000pan($idx4)), log10($rest3micronflux($idx4)), undef, $dlogrest3micronflux_err($idx4), {Symbol=>17,Colour=>'orange',SymbolSize=>2});

$idx1 = which( $rest3micronflux > $rest3micronflux_err & $spflag == 0 & $sfrOIIpan > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);
$idx2 = which( $rest3micronflux > $rest3micronflux_err & $spflag == 1 & $sfrOIIpan > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);
$idx3 = which( $rest3micronflux > $rest3micronflux_err & $spflag == 2 & $sfrOIIpan > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);
$idx4 = which( $rest3micronflux > $rest3micronflux_err & $spflag == -1 & $sfrOIIpan > 0 & ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);

$win->errb(log10($sfrOIIpan($idx1)), log10($rest3micronflux($idx1)), undef, $dlogrest3micronflux_err($idx1), {Symbol=>18,Colour=>'blue',SymbolSize=>2});
$win->errb(log10($sfrOIIpan($idx2)), log10($rest3micronflux($idx2)), undef, $dlogrest3micronflux_err($idx2), {Symbol=>18,Colour=>'green',SymbolSize=>2});
$win->errb(log10($sfrOIIpan($idx3)), log10($rest3micronflux($idx3)), undef, $dlogrest3micronflux_err($idx3), {Symbol=>18,Colour=>'red',SymbolSize=>2});
$win->errb(log10($sfrOIIpan($idx4)), log10($rest3micronflux($idx4)), undef, $dlogrest3micronflux_err($idx4), {Symbol=>18,Colour=>'orange',SymbolSize=>2});

$UL3micronflux = log10($rest3micronflux_err) + $dlogrest3micronflux_err;

#plot fit

$sfr = 10**(pdl(-4,-3,-2,1,0,1,2,3)); 
$r3mu = $int * $sfr**($slope);
$win->line(log10($sfr),log10($r3mu));

#use order of magnitude approximation Lexcess in Lband = 0.642 * (0.15) * sfr * lifetime
$Lexcess = 350; #calculated in circumstellar/disksize.pl
$flux3micron = $Lsun / 3e4 * $sfr * 1e6 * $Lexcess; #in units of W/A
$win->line(log10($sfr),log10($flux3micron),{Color=>green,Linestyle=>solid});
$win->legend(["Red Evolved","Intermediate","Star Forming","Mixed"],-1.0,32.5,{Colour=>['Red','Green','Blue','orange'],Symbol=>[17,17,17,17],Charsize=>1.1,SymbolSize=>1.5,TextShift=>-1.65});
$win->text("Toy CS Model",-0.9,29.6,{Angle=>30,Color=>green,Charsize=>1});

$win->close();


$opt = {Device=> "$ENV{HOME}/Software/paper/figures/AV_vs_mass.ps/vcps"};
$win = PDL::Graphics::PGPLOT::Window->new($opt);
$win->env(9.7,11.8,-0.1,2.1,{XTitle=>'Log M (M\d\(2281)\u)',YTitle=>'A\dv\u',AspectRatio=>1});


$idx1 = which( $z > 1.2 & $spflag == 0 &  ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);
$idx2 = which( $z > 1.2 & $spflag == 1 &  ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);
$idx3 = which( $z > 1.2 & $spflag == 2 &  ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);
$idx4 = which( $z > 1.2 & $spflag == -1 &  ($Conf > 1 | $z<9) & $Conf < 7 & $blendflag==0 & $IRACflag == 0);

$win->errb($logM($idx1), $AV($idx1), $logMerr($idx1), $AVerr($idx1), {Symbol=>17,Colour=>'blue',SymbolSize=>2});
$win->errb($logM($idx2), $AV($idx2), $logMerr($idx2), $AVerr($idx2), {Symbol=>17,Colour=>'green',SymbolSize=>2});
$win->errb($logM($idx3), $AV($idx3), $logMerr($idx3), $AVerr($idx3), {Symbol=>17,Colour=>'red',SymbolSize=>2});
$win->errb($logM($idx4), $AV($idx4), $logMerr($idx4), $AVerr($idx4), {Symbol=>17,Colour=>'orange',SymbolSize=>2});

#add line from Pannella et al

$logMpan = sequence(5) + 8.5;

#Pannella et al use a salpeter IMF
 
$AVpan = (peidust('SMC',5500) / peidust('SMC',1500)) * (4.07 * ($logMpan + log10(0.55) ) - 39.32);

$win->line($logMpan,$AVpan,{LineStyle=>'dashed',Color=>'blue'});
$win->close();
