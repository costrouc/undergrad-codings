// Fragment Shader

varying vec3 ec_vposition;
varying vec3 ts_veye;
varying vec2 ts_uvcoord;
varying mat3 TBN;

uniform sampler2D colormap;
uniform sampler2D normalmap;
uniform sampler2D bumpmap;

//Was planning to add enviroment maps...
uniform sampler2D diffuse_irr_map;
uniform sampler2D specular_irr_map;

void main()
{
  int numLights = 3;
  vec3 P,N,L,V,H;
  vec4 diffuse_color = gl_FrontMaterial.diffuse;
  vec4 specular_color = gl_FrontMaterial.specular;
  vec4 ambient_color = gl_FrontMaterial.ambient;
  float shininess = gl_FrontMaterial.shininess;
  
  //Color map
  vec3 texColor = texture2D(colormap, ts_uvcoord).rgb;
  
  //Normal map
  vec3 ts_vnormal = texture2D(normalmap, ts_uvcoord).rgb;
  ts_vnormal.rg = 2.0 * ts_vnormal.rg - 1.0;
  
  P = ec_vposition;
  N = normalize(ts_vnormal);
  V = normalize(ts_veye);
  
  gl_FragColor = vec4(0.0, 0.0, 0.0, 0.0);

  int i;
  for(i = 0; i < numLights; i++) {		
    L = TBN * normalize(gl_LightSource[i].position.xyz - P);
    H = normalize(L+V);

    diffuse_color = 0.0 * diffuse_color + 1.0 * vec4(texColor,1.0);
    diffuse_color *= max(dot(N,L),0.0);
    specular_color *= pow(max(dot(H,N),0.0),shininess);
	
    gl_FragColor += ambient_color*gl_LightSource[i].ambient;
    gl_FragColor += diffuse_color*gl_LightSource[i].diffuse;
    gl_FragColor += specular_color*gl_LightSource[i].specular;
  }
}
