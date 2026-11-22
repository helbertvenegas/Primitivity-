# Verificaciones del artículo: versión completada 0.3.0

Resultado: **2.900 de 2.900 comprobaciones aprobadas**, sin fallos.

El índice recorre las 108 etiquetas del manuscrito y exige certificados para 57 bloques de cálculo. No faltan bloques obligatorios.

El verificador independiente aprobó **1.840 comprobaciones**, utilizando únicamente la biblioteca estándar de Python.

El cuaderno ejecutó sus **31 celdas de código** mediante nbclient y el kernel oficial de IPython en modo local (`inprocess`), en un directorio vacío y sin importar el repositorio. Se conservaron las salidas de Jupyter y no hubo errores.

El entorno no permite abrir los canales de comunicación de un kernel externo. Ese modo queda disponible en el ejecutor y configurado en GitHub Actions, pero no se presenta como probado aquí. Colab solicitó iniciar sesión en Google y su ejecución remota queda pendiente.

## Pendientes computacionales cerrados

- Inducción general de pesos: proyectores simbólicos, obstrucción mediante Omega y recubrimiento de índices para todo n.
- Fórmulas generales de grado PI: separación de toros, todas las clases de paridad, parametrizaciones centrales completas y equivalencias de ideales enteros para D/G/K/L.
- Segundo verificador: aritmética racional propia, otra forma normal de palabras, testigos unimodulares de Smith y de núcleos completos; sin importar SymPy ni el paquete principal.
- Trazabilidad por afirmación: hipótesis, identidad comprobada, identificadores de certificados y dependencias teóricas.

## Cobertura del paquete

- Los modelos de las 21 familias con fórmula cerrada de grado PI, con matrices de conmutación y de acción obtenidas de las relaciones.
- El cálculo exacto de la red central, el grupo de acciones sobre el centro y los factores del producto cruzado. Se admiten órdenes arbitrarios, sin el límite de las muestras anteriores.
- Las identidades de Theta y Z1,Z2 para F,I,Q,R,V; las anti-isomorfías y la separación simbólica de los toros en factores.
- Elementos centrales racionales y por órbitas, normalidad de localizadores, obstrucciones a acciones interiores, Omega/Xi y matrices de pesos.
- Certificados generales de Dickson, recurrencias y una reducción de módulos de rango nueve para C.
- Una correspondencia entre comprobaciones y lemas/ecuaciones, con las dependencias teóricas identificadas.

## Corrección detectada en el manuscrito

En la matriz de L, en el orden w,y1,y1^-1*y2,u, la primera entrada superior debe ser **-(2 eta+phi)**. La relación y1*w=h^2*f*w*y1 exige ese signo con la convención usada.
El cambio conserva el Pfaffiano y los invariantes alternantes de Smith, por lo que **no cambia la fórmula del grado PI**. El cuaderno contiene una comprobación que detecta el signo antiguo. No se modificó silenciosamente el manuscrito.

## Alcance exacto

Los certificados simbólicos se distinguen de las pruebas de regresión con órdenes finitos. La decisión PI aplica la tabla del manuscrito; el cálculo del grado se obtiene del toro y sus acciones, y después se compara con la fórmula.
Las pruebas abstractas sobre anillos de división, independencia, ideales y Nullstellensatz siguen siendo argumentos matemáticos del artículo. No se afirma una formalización en un asistente de pruebas.
El grado exacto de C y los casos de primitividad que el artículo deja abiertos se conservan como preguntas abiertas.

## Archivos para revisar

- `docs/MATRICES.md`: todas las matrices de definición.
- `docs/TORUS_ACTION_MATRICES.md`: matrices de los 21 modelos de grado.
- `reports/coverage.md` y `coverage.json`: correspondencia completa con las etiquetas.
- `docs/SYMBOLIC_CERTIFICATES.md`: enlace entre certificados y argumentos generales.
- `docs/UNBOUNDED_PROOFS.md`: inducción y demostraciones generales para todos los órdenes.
- `independent/check.py`: verificador independiente; `reports/certificates.json`: datos y testigos.
- `reports/notebook_execution.json`: constancia del modo de ejecución y las huellas de los módulos.
- `docs/MANUSCRIPT_CORRECTIONS.md`: corrección de L lista para incorporar.

Entorno: Python 3.12.14; SymPy 1.14.0.
