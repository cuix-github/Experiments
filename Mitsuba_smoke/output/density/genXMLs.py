import xml.etree.ElementTree as etree
import xml.dom.minidom
import subprocess
import numpy as np

total_frames = 240;
for i in range(0, total_frames):
    print "Generating %d frame scene XML file" % (i)

    scene = etree.Element("scene", version="0.5.0")
    integrator = etree.SubElement(scene, "integrator", type="volpath_simple")
    etree.SubElement(integrator, "integer", name="maxDepth", value="8")
    
    # Medium
    medium = etree.SubElement(scene, "medium", type="heterogeneous", id="smoke")
    method = etree.SubElement(medium, "string", name="method", value="woodcock")
    density = etree.SubElement(medium, "volume", name="density", type="gridvolume")
    etree.SubElement(density, "string", name="filename", value="density-%04d.vol" % i)
    albedo = etree.SubElement(medium, "volume", name="albedo", type="constvolume")
    etree.SubElement(albedo, "spectrum", name="value", value="0.9")
    etree.SubElement(medium, "float", name="scale", value="100")
    
    # Interior boundary
    bounds = etree.SubElement(scene, "shape", type="obj")
    etree.SubElement(bounds, "string", name="filename", value="../bounds.obj")
    etree.SubElement(bounds, "ref", name="interior", id="smoke")
    
    # Ground
    floor = etree.SubElement(scene, "shape", type="obj")
    etree.SubElement(floor, "string", name="filename", value="../plane.obj")
    material = etree.SubElement(floor, "bsdf", type="diffuse")
    etree.SubElement(material, "rgb", name="diffuseReflectance", value=".2, .2, .3")
    floor_transform = etree.SubElement(floor, "transform", name="toWorld")
    etree.SubElement(floor_transform, "translate", y=".48")
    
    # Sensor (Camera settings)
    sensor = etree.SubElement(scene, "sensor", type="perspective")
    etree.SubElement(sensor, "float", name="focusDistance", value="1.25668")
    etree.SubElement(sensor, "float", name="fov", value="45")
    etree.SubElement(sensor, "string", name="fovAxis", value="x")
    sensor_transform = etree.SubElement(sensor, "transform", name="toWorld")
    etree.SubElement(sensor_transform, "scale", x="-1")
    etree.SubElement(sensor_transform, "lookat", target="-0.166029, 0.148984, -0.537402",
                                                 origin="-0.61423, 0.154197, -1.43132",
                                                 up="-0.000640925, -0.999985, -0.0055102")
    sensor_sampler = etree.SubElement(sensor, "sampler", type="ldsampler")
    etree.SubElement(sensor_sampler, "integer", name="sampleCount", value="512")
    film_settings = etree.SubElement(sensor, "film", type="hdrfilm")
    etree.SubElement(film_settings, "integer", name="height", value="576")
    etree.SubElement(film_settings, "integer", name="width", value="768")
    etree.SubElement(film_settings, "rfilter", type="gaussian")

    # Illumination
    light = etree.SubElement(scene, "shape", type="sphere")
    etree.SubElement(light, "point", name="center", x="0", y="-2", z="-1")
    etree.SubElement(light, "float", name="radius", value=".2")
    emitter_settings = etree.SubElement(light, "emitter", type="area")
    etree.SubElement(emitter_settings, "spectrum", name="radiance", value="400")

    rough_string = etree.tostring(scene, "utf-8")
    reparsed = xml.dom.minidom.parseString(rough_string)

    reparsed_pretty = reparsed.toprettyxml(indent=" " * 4)

    with open("frame_%04d.xml" % i, "w") as frame_xml:
        frame_xml.write(reparsed_pretty)