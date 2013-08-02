// Takes in account 3 lights
// Phong lighting in eye coordinates.

// These are set by the .vert code, interpolated. 
varying vec3 ec_vnormal;
varying vec4 ec_vposition;

void main()
{
  vec3 N,L,V,H;
  vec4 P;
  float shininess = gl_FrontMaterial.shininess;
  vec4 diffuse_color;
  vec4 specular_color;
  vec4 ambient_color = gl_FrontMaterial.ambient;

  int i; 
  int numLights = 3;

  P = normalize(ec_vposition);
  N = normalize(ec_vnormal);

  V = normalize(P.xyz); //eye position is (0,0,0)!;

  gl_FragColor = vec4(0.0,0.0,0.0,0.0);

  for(i = 0; i < numLights; i++) {
	
	diffuse_color = gl_FrontMaterial.diffuse;
	specular_color = gl_FrontMaterial.specular;
	
	L = normalize(gl_LightSource[i].position - P).xyz;
	H = normalize(L+V);
	
	diffuse_color *= max(dot(N,L),0.0);
    specular_color *= pow(max(dot(H,N),0.0),shininess);
	
	gl_FragColor += ambient_color*gl_LightSource[i].ambient;
	gl_FragColor += diffuse_color*gl_LightSource[i].diffuse;
	gl_FragColor += specular_color*gl_LightSource[i].specular;
  }
}