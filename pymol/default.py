# BallnStick creates a ball and stick representation of an object 
# Add_VDW creates a copy of an object with full-sized, transparent spheres
# Bondi VDW values added below to override default Pymol settings

from pymol import cmd
from pymol import util
from pymol import stored
import re

#Color Maps
viridis = [
    [0.992, 0.905, 0.145],
    [0.368, 0.788, 0.384],
    [0.129, 0.568, 0.549],
    [0.231, 0.321, 0.545],
    [0.266, 0.003, 0.329]
]


# Bondi VDW values 
cmd.alter("elem Ac", "vdw=2.00")
cmd.alter("elem Al", "vdw=2.00")
cmd.alter("elem Am", "vdw=2.00")
cmd.alter("elem Sb", "vdw=2.00")
cmd.alter("elem Ar", "vdw=1.88")
cmd.alter("elem As", "vdw=1.85")
cmd.alter("elem At", "vdw=2.00")
cmd.alter("elem Ba", "vdw=2.00")
cmd.alter("elem Bk", "vdw=2.00")
cmd.alter("elem Be", "vdw=2.00")
cmd.alter("elem Bi", "vdw=2.00")
cmd.alter("elem Bh", "vdw=2.00")
cmd.alter("elem B ", "vdw=2.00")
cmd.alter("elem Br", "vdw=1.85")
cmd.alter("elem Cd", "vdw=1.58")
cmd.alter("elem Cs", "vdw=2.00")
cmd.alter("elem Ca", "vdw=2.00")
cmd.alter("elem Cf", "vdw=2.00")
cmd.alter("elem C ", "vdw=1.70")
cmd.alter("elem Ce", "vdw=2.00")
cmd.alter("elem Cl", "vdw=1.75")
cmd.alter("elem Cr", "vdw=2.00")
cmd.alter("elem Co", "vdw=2.00")
cmd.alter("elem Cu", "vdw=1.40")
cmd.alter("elem Cm", "vdw=2.00")
cmd.alter("elem Ds", "vdw=2.00")
cmd.alter("elem Db", "vdw=2.00")
cmd.alter("elem Dy", "vdw=2.00")
cmd.alter("elem Es", "vdw=2.00")
cmd.alter("elem Er", "vdw=2.00")
cmd.alter("elem Eu", "vdw=2.00")
cmd.alter("elem Fm", "vdw=2.00")
cmd.alter("elem F ", "vdw=1.47")
cmd.alter("elem Fr", "vdw=2.00")
cmd.alter("elem Gd", "vdw=2.00")
cmd.alter("elem Ga", "vdw=1.87")
cmd.alter("elem Ge", "vdw=2.00")
cmd.alter("elem Au", "vdw=1.66")
cmd.alter("elem Hf", "vdw=2.00")
cmd.alter("elem Hs", "vdw=2.00")
cmd.alter("elem He", "vdw=1.40")
cmd.alter("elem Ho", "vdw=2.00")
cmd.alter("elem In", "vdw=1.93")
cmd.alter("elem I ", "vdw=1.98")
cmd.alter("elem Ir", "vdw=2.00")
cmd.alter("elem Fe", "vdw=2.00")
cmd.alter("elem Kr", "vdw=2.02")
cmd.alter("elem La", "vdw=2.00")
cmd.alter("elem Lr", "vdw=2.00")
cmd.alter("elem Pb", "vdw=2.02")
cmd.alter("elem Li", "vdw=1.82")
cmd.alter("elem Lu", "vdw=2.00")
cmd.alter("elem Mg", "vdw=1.73")
cmd.alter("elem Mn", "vdw=2.00")
cmd.alter("elem Mt", "vdw=2.00")
cmd.alter("elem Md", "vdw=2.00")
cmd.alter("elem Hg", "vdw=1.55")
cmd.alter("elem Mo", "vdw=2.00")
cmd.alter("elem Nd", "vdw=2.00")
cmd.alter("elem Ne", "vdw=1.54")
cmd.alter("elem Np", "vdw=2.00")
cmd.alter("elem Ni", "vdw=1.63")
cmd.alter("elem Nb", "vdw=2.00")
cmd.alter("elem N ", "vdw=1.55")
cmd.alter("elem No", "vdw=2.00")
cmd.alter("elem Os", "vdw=2.00")
cmd.alter("elem O ", "vdw=1.52")
cmd.alter("elem Pd", "vdw=1.63")
cmd.alter("elem P ", "vdw=1.80")
cmd.alter("elem Pt", "vdw=1.72")
cmd.alter("elem Pu", "vdw=2.00")
cmd.alter("elem Po", "vdw=2.00")
cmd.alter("elem K ", "vdw=2.75")
cmd.alter("elem Pr", "vdw=2.00")
cmd.alter("elem Pm", "vdw=2.00")
cmd.alter("elem Pa", "vdw=2.00")
cmd.alter("elem Ra", "vdw=2.00")
cmd.alter("elem Rn", "vdw=2.00")
cmd.alter("elem Re", "vdw=2.00")
cmd.alter("elem Rh", "vdw=2.00")
cmd.alter("elem Rb", "vdw=2.00")
cmd.alter("elem Ru", "vdw=2.00")
cmd.alter("elem Rf", "vdw=2.00")
cmd.alter("elem Sm", "vdw=2.00")
cmd.alter("elem Sc", "vdw=2.00")
cmd.alter("elem Sg", "vdw=2.00")
cmd.alter("elem Se", "vdw=1.90")
cmd.alter("elem Si", "vdw=2.10")
cmd.alter("elem Ag", "vdw=1.72")
cmd.alter("elem Na", "vdw=2.27")
cmd.alter("elem Sr", "vdw=2.00")
cmd.alter("elem S ", "vdw=1.80")
cmd.alter("elem Ta", "vdw=2.00")
cmd.alter("elem Tc", "vdw=2.00")
cmd.alter("elem Te", "vdw=2.06")
cmd.alter("elem Tb", "vdw=2.00")
cmd.alter("elem Tl", "vdw=1.96")
cmd.alter("elem Th", "vdw=2.00")
cmd.alter("elem Tm", "vdw=2.00")
cmd.alter("elem Sn", "vdw=2.17")
cmd.alter("elem Ti", "vdw=2.00")
cmd.alter("elem W ", "vdw=2.00")
cmd.alter("elem U ", "vdw=1.86")
cmd.alter("elem V ", "vdw=2.00")
cmd.alter("elem Xe", "vdw=2.16")
cmd.alter("elem Yb", "vdw=2.00")
cmd.alter("elem Y ", "vdw=2.00")
cmd.alter("elem Zn", "vdw=1.39")
cmd.alter("elem Zr", "vdw=2.00")
cmd.rebuild()

# workspace settings
cmd.set("internal_gui_width", 400)
cmd.bg_color("white")

cmd.set("ray_trace_mode", 1)
cmd.set("ray_texture", 0)
cmd.set("ray_opaque_background", "off")

cmd.set('dash_color', 'gray')
cmd.set("dash_gap",0.15)
cmd.set("dash_radius",0.035)

cmd.space("cmyk")

cmd.set('label_color', 'yellow')
cmd.set('label_outline_color', 'black')
cmd.set('label_font_id', 7)
cmd.set('label_size', 30)
###################
def default_env(arg1):
	cmd.set("ambient", 0.4) #amount of ambient light
	cmd.set("shininess", 25) #how much the object shines <the greater it becames opaque>
	cmd.set("reflect", 0.05) #amount of light reflection
	cmd.set("orthoscopic", 0)
	cmd.set("transparency", 0.5)
	cmd.set("antialias", 3)
	cmd.set("spec_count", 5)
	cmd.set("specular", 1.5)

def default_colors(arg1):
    util.cbaw(arg1)
    cmd.color("gray35", f"elem C and {arg1}")
    cmd.color("gray95", f"elem H and {arg1}")
    cmd.color("br9", f"elem 0 and {arg1}")
    cmd.color("br0", f"elem N and {arg1}")

    cmd.set('stick_transparency', 0, arg1)
    cmd.set('stick_color', 'default', arg1)
    cmd.set('sphere_transparency', 0, arg1)

def ball_and_stick(arg1='all'):
    default_colors(arg1)
    default_env(arg1)
    cmd.show("sticks", arg1)
    cmd.show("spheres", arg1)
    cmd.hide("nonbonded", arg1)
    cmd.hide("lines", arg1)
    cmd.hide("labels")

    cmd.set("stick_radius",0.2, arg1)
    cmd.set("stick_h_scale",1.0, arg1) 
    cmd.set("sphere_scale",0.25, arg1)
    cmd.set("sphere_scale",0.25, 'elem H')
cmd.extend("bns", ball_and_stick)

def houkmol(arg1='all'):
    ball_and_stick(arg1)
    default_env(arg1)
    
    cmd.color("gray70", f"elem C and {arg1}")
    cmd.color("red", f"elem O and {arg1}")
    cmd.color("tv_blue", f"elem N and {arg1}")
    
    cmd.set("stick_color",'black', arg1)
    cmd.set("stick_radius",0.1, arg1)
    cmd.set("stick_h_scale",1.0, arg1) 
    
    cmd.set('sphere_scale',0.20, arg1)
    cmd.set('sphere_scale',0.15, f'elem H and {arg1}')
    
    cmd.set("shininess", 100)
    cmd.set("ambient", 0.6)
    
cmd.extend("hmol", houkmol)

def quick_overlay(entry, color='blue'):
    cmd.color(color,entry)
    cmd.set('stick_transparency', 0.5, entry)
    cmd.set('sphere_transparency', 0.5, entry)
cmd.extend('quickover', quick_overlay)


def add_VDW(arg1):
    cmd.copy(arg1+"_vdw", arg1)
    cmd.set("sphere_scale",1.0, arg1+"_vdw and elem H")
    cmd.rebuild()
    cmd.set("sphere_scale", 1, arg1+"_vdw")
    cmd.hide("nonbonded", arg1+"_vdw")
    cmd.hide("lines", arg1+"_vdw")
    cmd.hide("sticks", arg1+"_vdw")
    cmd.set("sphere_transparency", 0.5, arg1+"_vdw")
cmd.extend("add_vdw", add_VDW)


def plot_cube(isovalue=0.004):
    cmd.set("internal_gui_width", 525)
    obj_list = cmd.get_names('objects')
    isovalue = float(isovalue)
    for cube in obj_list:
        orbitalName = 'orb-'+cube
        positiveOrbital = cube+'+'
        negativeOrbital = cube+'-'
        
        cmd.isomesh(positiveOrbital,cube,isovalue*1)
        cmd.color("blue",positiveOrbital) 
        cmd.isomesh(negativeOrbital,cube,isovalue*-1)
        cmd.color("red",negativeOrbital)

        for orbital in (positiveOrbital,negativeOrbital):
            cmd.group(orbitalName, orbital)
cmd.extend("plot_cube", plot_cube)


def bond_between(element1, element2, cutoff=2.3):
    visible_entries = cmd.get_object_list(selection='visible')
 
    for entry in visible_entries:

        bonds = cmd.find_pairs(
            f"element {element1} and {entry}", 
            f"element {element2} and {entry}", 
            cutoff=cutoff)
        print(f"{entry} has {len(bonds)} bonds between {element1} and {element2}")
        
        for bond in bonds:
            entry = bond[0][0]
            element1_id = bond[0][1]
            element2_id = bond[1][1]
            cmd.bond(f"id {element1_id} and {entry}", f"id {element2_id} and {entry}")
cmd.extend("bondbetween", bond_between)            

def nci(arg1, isovalue=0.3):
	# nci.py, a tiny script to display plots from Nciplot in PyMOL
    densf = arg1+"-dens"
    gradf = arg1+"-grad"
    cmd.isosurface("grad",gradf, isovalue)
    cmd.ramp_new("ramp", densf, [-5,5], 'rainbow')
    cmd.set("surface_color", "ramp", "grad")
    cmd.set('transparency', 0, 'grad')
    cmd.set('two_sided_lighting',value=1)
cmd.extend( "nci", nci );

def electrostatic_potential_surface(arg1, isovalue=0.04, scale=0.5):
    densf = arg1+"_dens"
    gradf = arg1+"_esp"
    min_surface, max_surface = [abs(float(scale))*-1, abs(float(scale))]
    cmd.isosurface("dens",densf, isovalue)
    cmd.ramp_new("ramp", gradf, [min_surface, max_surface], 'rainbow')
    cmd.set("surface_color", "ramp", "dens")
    cmd.set('transparency', 0.50, 'dens')
    cmd.set('two_sided_lighting',value=1)
cmd.extend( "elpot", electrostatic_potential_surface );

def group_visible(groupname='test', include_measurements=False):
    visible_entries = cmd.get_object_list(selection='visible')
    if include_measurements:
        all_objects = cmd.get_names(type='all')
        for item in all_objects:
            dist_obj = re.search('dist([0-9]{1,4})', item)
            if dist_obj:
                dist_obj_id = dist_obj.group(1)
                dist_obj_newname = f"{visible_entries[0]}_dista{dist_obj_id}"
                cmd.set_name(item, dist_obj_newname)
                visible_entries.append(dist_obj_newname)
    cmd.group(groupname, ' '.join(visible_entries))
cmd.extend("group_visible", group_visible);

def show_buried_volume(entry_name, sphere_radius=3.5):
    cmd.set_name(entry_name, "bv_atoms")
    cmd.copy(entry_name, "bv_atoms")
    cmd.alter("sele", f"vdw={sphere_radius}")
    cmd.set("sphere_scale", 1, "sele")
    cmd.set("sphere_transparency", 0.50, "sele")
cmd.extend("bv", show_buried_volume);

def print_stuff():
    visible_entries = cmd.get_names(type='all', enabled_only=1)
    print(visible_entries)
cmd.extend("pstuff", print_stuff);

def align_visible(reference_entry_id=0):
    visible_entries = cmd.get_object_list(selection='visible')
    reference_entry = visible_entries[reference_entry_id]
    mobile_entries = [i for i in visible_entries if i != reference_entry]
    for entry in mobile_entries:
        cmd.align(entry, reference_entry, cycles=200)
cmd.extend("align_visible", align_visible);

def save_img(filename='default', ray_mode=1):
    '''
    Export the current workspace as PNG using Ray-Tracing. If filename is not present
    the image is named default.png, except if only one element is visible in the workspace.
    In this case the filename is the same as the name of the visible entry.
    '''
    visible_entries = cmd.get_object_list(selection='visible')
    if filename == 'default' and len(visible_entries) == 1:
        filename = visible_entries[0]

    screen_size = cmd.get_viewport(output=1, quiet=1)
    canvas_width = screen_size[0]
    canvas_height = screen_size[1]

    cmd.png(filename, canvas_width*3, canvas_height*3, dpi=300, ray=ray_mode)
cmd.extend("expimg", save_img)



