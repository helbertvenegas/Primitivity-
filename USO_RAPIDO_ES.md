# Uso rápido y alcance de la versión completada

El paquete acompaña el artículo de James Gómez y Helbert Venegas sobre las 26
familias de extensiones dobles de Ore de tipo (14641).

## Ejecutar en Colab

1. Abre Google Colab y elige Archivo > Subir cuaderno.
2. Selecciona `Double_Ore_14641_PI_verification.ipynb` y ejecuta todas las celdas.
3. Cambia `family_name` para inspeccionar una familia A--Z, o cambia la etiqueta
   del lema para consultar sus certificados.

El cuaderno contiene su código; no requiere Drive, un repositorio externo ni GPU.
Sus explicaciones principales están en inglés para acompañar el artículo.

## Validación de esta versión

La versión 0.3.0 añade certificados para cualquier grado u orden, un segundo
verificador sin SymPy, hipótesis por afirmación y ejecución mediante nbclient
y el kernel oficial de IPython en modo local (`inprocess`). El informe de
ejecución conserva los conteos, el tipo de kernel y las huellas del código.
Colab requiere iniciar sesión en Google; esa ejecución remota queda pendiente.
El modo de kernel externo está disponible en el ejecutor y en GitHub Actions.

## Qué se completó

Además de las matrices de definición, se incluyen los modelos de toros y las
acciones para las 21 familias con fórmula cerrada de grado PI. Las matrices se
obtienen de las relaciones mediante cálculo exacto con inversos. El programa
calcula la red central y el grupo de acciones sobre el centro y compara el
resultado con la fórmula del manuscrito, sin un límite fijo para el orden.

También se incluyeron las identidades de Theta y de Z1,Z2, las anti-isomorfías,
los elementos normales de localización, certificados centrales racionales y
por órbitas, la paridad que obstruye acciones interiores, Omega/Xi, matrices de
pesos, los certificados generales de Dickson y la base de rango nueve para C.

El registro recorre las 108 etiquetas del manuscrito. Distingue cálculos
comprobados, dependencias teóricas, pruebas finitas adicionales y preguntas
abiertas. Las demostraciones abstractas de teoría de ideales y anillos de
división siguen siendo argumentos matemáticos del artículo; no se anuncia una
formalización en un asistente de pruebas.

## Corrección encontrada en el manuscrito

En la matriz de conmutación de L, ordenada como w,y1,y1^-1*y2,u, debe escribirse
-(2*eta+phi) en la primera entrada superior. La relación y1*w=h^2*f*w*y1 lo
comprueba. El cambio conserva los invariantes de Smith y no cambia la fórmula
del grado PI. La nota `MANUSCRIPT_CORRECTIONS.md` contiene el reemplazo exacto.

## Ejemplo con parámetros exactos

```python
from sympy import Rational
from ore14641 import f,h
from ore14641.degrees import degree_from_parameters

resultado = degree_from_parameters('Z',{f:3,h:-Rational(1,2)})
print(resultado['algebra_pi_degree'])  # 6
```

Para C se devuelven las tres posibilidades de su grado, porque su determinación
exacta sigue abierta en el manuscrito. Las preguntas restantes de primitividad
tampoco se presentan como resueltas por el código.

## Archivos para revisar

- `docs/MATRICES.md`: matrices Sigma,M y det(sigma) de A--Z.
- `docs/TORUS_ACTION_MATRICES.md`: matrices de conmutación y de acción de los modelos.
- `docs/SYMBOLIC_CERTIFICATES.md`: enlace entre los certificados y los argumentos generales.
- `docs/UNBOUNDED_PROOFS.md`: argumentos generales de pesos, paridad y grados.
- `reports/certificates.json`: datos para el verificador independiente.
- `reports/independent_verification.json`: resultado de la segunda comprobación.
- `reports/coverage.json`: hipótesis y certificados por afirmación del artículo.
- `reports/verification.json`: todos los residuos y certificados.

Para GitHub, sube el contenido del ZIP y conserva `.github/workflows/verify.yml`.
Nombre sugerido: `double-ore-14641-pi-verification`. La publicación y la elección
de una licencia corresponden a una etapa posterior; este paquete no incluye
el PDF editorial ni el manuscrito completo.
