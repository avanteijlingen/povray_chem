global_settings {
	ambient_light rgb <0.200000002980232, 0.200000002980232, 0.200000002980232>
	max_trace_level 15
}

background { color rgb <1,1,1> }

camera {
	perspective
	location <0, 0, 20>
	angle 40
	up <0, 1, 0>
	right <1, 0, 0> * 1
	direction <0, 0, -1> }

light_source {
	<31.7157010466102, 27.7512375296575, 39.6446257175118>
	color rgb <1, 1, 1>
	fade_distance 79.2892514350236
	fade_power 0
	parallel
	point_at <-31.7157010466102, -27.7512375296575, -39.6446257175118>
}

light_source {
	<-31.7157010466102, 27.7512375296575, -19.8223128587559>
	color rgb <0.300000011920929, 0.300000011920929, 0.300000011920929>
	fade_distance 79.2892514350236
	fade_power 0
	parallel
	point_at <31.7157010466102, -27.7512375296575, 19.8223128587559>
}

#default {
	finish {ambient .8 diffuse 1 specular 1 roughness .005 metallic 0.5}
}


