varying vec3 normal,lightDir;

void main()
{
    normal = normalize(gl_NormalMatrix * gl_Normal);
    
    vec3 worldPos = vec3(gl_ModelViewMatrix * gl_Vertex);
    lightDir = normalize(vec3(gl_LightSource[0].position) - worldPos);
    
    gl_Position = ftransform();
}
