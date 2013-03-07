// Phong lighting in eye coordinates.

// These are set by the .vert code, interpolated. 
varying vec3 ec_vnormal, ec_vposition;

uniform sampler2D textureMap;
uniform sampler2D normalMap;
uniform sampler2D dispMap;

void main()
{
  vec3 P,N,L,V,H;
  vec4 diffuse_color = gl_FrontMaterial.diffuse;
  vec4 specular_color = gl_FrontMaterial.specular;
  float shininess = gl_FrontMaterial.shininess;

  vec3 normal = 2.0 * texture2D(normalMap, gl_TexCoord[0].st).rgb -1.0;
  
  P = ec_vposition;
  //N = normalize(ec_vnormal);
  N = normalize(normal + ec_vnormal);
  L = normalize(gl_LightSource[0].position.xyz - P);
  V = normalize(P); //eye position is (0,0,0)!;
  H = normalize(L+V);

  vec3 texColor = texture2D(textureMap,gl_TexCoord[0].st).rgb;
  //Blend base diffuse color with texture color
  diffuse_color = 0.0 * diffuse_color + 1.0 * vec4(texColor,1.0);
  diffuse_color *= max(dot(N,L),0.0);
  specular_color *= pow(max(dot(H,N),0.0),shininess);
  gl_FragColor = diffuse_color + specular_color;
  //gl_FragColor = vec4(texColor,1.0);
}
