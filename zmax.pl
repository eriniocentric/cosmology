use lib "$ENV{HOME}/Software";
$ENV{DATADIR} = "$ENV{HOME}/Software/filters";

use PDL;
use PDL::NiceSlice;
use KGB::Cosmology;
use KGB::PegUtils;
use KGB::SpecUtils;
use KGB::Dust;

$ABflag = 0; #use Vega K mag to compare for zmax
$OmegaM = 0.3; $H0 = 70;
$Law = "SMC"; # Dust Law MW/SMC
$nebfac = 2.0; # Extra factor of extinction for neb lines
$pc = 3.086E16; # One parsec in meters

open IN, "$ENV{HOME}/Software/Masses/GDDSVIzK12/SEDparameters.dat" or die "Can't open SED file";
open OUT, ">$ENV{HOME}/Software/Masses/GDDSVIzK12/zmax.dat" or die "can't open zmax file"; 
print OUT "# id      zmax\n";
select OUT; $| = 1; select STDOUT; $| = 1;

while (<IN>) {
    next if /^#/;
    @v = split;
    $id = $v[0]; $zsp = $v[1]; $specfile = $v[2]; $itime = $v[3]; $fac = $v[4]; $AVV =$v[5];

    print "Getting SED from $specfile\n";

    ($t2,$dummy, $wav,$spec,$emspec) = read_peg_spec($specfile,{SPLIT_EMISSION=>1}); # Note #spec returned in W/A
    $spec *= 1E10; $emspec *= 1E10; # Convert to 1E10 mass galaxy to avoid rounding
         
    $spec2 = $spec->dice_axis(1,$itime); 
    $emspec2 = $emspec->dice_axis(1,$itime); 
           
    $attn = peidust($Law, $wav)/peidust($Law, 5500);
    $spec2 = $spec2 * 10**( -0.4 * $attn * $AVV ) + $emspec2 * 10**(-0.4 * $attn * $nebfac * $AVV);
    $z = 0;
    $K = 0;
    while ($K <= 20.6) {
	$z+= 0.001;
	($wz, $sz) = redshift_spectra($wav, $spec2, pdl($z));
	
	# Coerce 1D             
	die "Something wrong \n" if (nelem($sz) ne $sz->getdim(0)) or (nelem($wz) ne $wz->getdim(0)); # Should never happen! 
	$wz = $wz->clump(-1)->copy; $sz = $sz->clump(-1)->copy;
	
	$DL = lumdist($z);
	$sz /= 4*$Pi * ($DL*1E6*$pc)**2; 
	
	# $sz now matches that used to produce the colors for this model
	$sz *= $fac; # Normalization
	$K = mag($wz,$sz,'sys_circi_K');
    }
    $zmax = $z;
    $zmax = 0.0 if ($zmax < $zsp);
    printf OUT "%10s %5.2f\n", $id,$zmax;
}
close IN; close OUT;
