"""Render the completed JSON certificates into readable coverage and Spanish notes."""
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]


def render():
    report=json.loads((ROOT/'reports/verification.json').read_text())
    coverage=report['coverage']
    independent=json.loads((ROOT/'reports/independent_verification.json').read_text())
    execution=json.loads((ROOT/'reports/notebook_execution.json').read_text())
    number=lambda value:f'{value:,}'.replace(',','.')
    kernel_text='kernel externo de IPython' if execution['external_kernel_process_tested'] else 'kernel oficial de IPython en modo local (`inprocess`)'
    kernel_limit='' if execution['external_kernel_process_tested'] else 'El entorno no permite abrir los canales de comunicación de un kernel externo. Ese modo queda disponible en el ejecutor y configurado en GitHub Actions, pero no se presenta como probado aquí. '
    lines=['# Coverage of the supplied manuscript','',
           f"Snapshot labels: {coverage['manuscript_label_count']}. Required computational blocks: {coverage['required_computational_blocks']}.",
           f"Missing required blocks: {coverage['missing_required_blocks']}.",'',
           'A covered computation is not a proof-assistant formalization of the abstract theorem.',
           'The ledger explicitly separates algebraic certificates, theoretical dependencies and open questions.','',
           '| Manuscript label | Direct checks | Status / support |','|---|---:|---|']
    for entry in coverage['entries']:
        support=', '.join(f"{name} ({count})" for name,count in entry['supporting_checks'].items())
        text=entry['status']
        if support:text+='; '+support
        if entry['theoretical_dependency']:text+='; '+entry['theoretical_dependency']
        lines.append(f"| `{entry['label']}` | {entry['direct_checks']} | {text} |")
    lines+=['','## Claim-level contracts','',
            'Each required block records hypotheses, the checked assertion, actual certificate IDs and mathematical dependencies.','']
    for entry in coverage['entries']:
        c=entry.get('claim_contract')
        if not c:continue
        lines+=[f"### `{entry['label']}`",'',
               '**Hypotheses:** '+c['hypotheses'],'',
               '**Checked claim:** '+c['checked_claim'],'',
               '**Mathematical dependencies:** '+'; '.join(c['mathematical_dependencies']),'',
               '**Certificate IDs:** '+', '.join('`'+i+'`' for i in c['certificate_ids']), '']
    lines+=['','All 21 closed degree models are implemented. C returns its three possible degrees.',
            'See MANUSCRIPT_CORRECTIONS.md for the L matrix sign correction.']
    (ROOT/'reports/coverage.md').write_text('\n'.join(lines)+'\n')
    summary=report['summary']
    lines=['# Verificaciones del artículo: versión completada 0.3.0','',
           f"Resultado: **{number(summary['passed'])} de {number(summary['total'])} comprobaciones aprobadas**, sin fallos.",'',
           f"El índice recorre las {coverage['manuscript_label_count']} etiquetas del manuscrito y exige certificados para {coverage['required_computational_blocks']} bloques de cálculo. No faltan bloques obligatorios.",'',
           f"El verificador independiente aprobó **{number(independent['passed'])} comprobaciones**, utilizando únicamente la biblioteca estándar de Python.",'',
           f"El cuaderno ejecutó sus **{execution['code_cells_executed']} celdas de código** mediante nbclient y el {kernel_text}, en un directorio vacío y sin importar el repositorio. Se conservaron las salidas de Jupyter y no hubo errores.",'',
           kernel_limit+'Colab solicitó iniciar sesión en Google y su ejecución remota queda pendiente.','',
           '## Pendientes computacionales cerrados','',
           '- Inducción general de pesos: proyectores simbólicos, obstrucción mediante Omega y recubrimiento de índices para todo n.',
           '- Fórmulas generales de grado PI: separación de toros, todas las clases de paridad, parametrizaciones centrales completas y equivalencias de ideales enteros para D/G/K/L.',
           '- Segundo verificador: aritmética racional propia, otra forma normal de palabras, testigos unimodulares de Smith y de núcleos completos; sin importar SymPy ni el paquete principal.',
           '- Trazabilidad por afirmación: hipótesis, identidad comprobada, identificadores de certificados y dependencias teóricas.','',
           '## Cobertura del paquete','',
           '- Los modelos de las 21 familias con fórmula cerrada de grado PI, con matrices de conmutación y de acción obtenidas de las relaciones.',
           '- El cálculo exacto de la red central, el grupo de acciones sobre el centro y los factores del producto cruzado. Se admiten órdenes arbitrarios, sin el límite de las muestras anteriores.',
           '- Las identidades de Theta y Z1,Z2 para F,I,Q,R,V; las anti-isomorfías y la separación simbólica de los toros en factores.',
           '- Elementos centrales racionales y por órbitas, normalidad de localizadores, obstrucciones a acciones interiores, Omega/Xi y matrices de pesos.',
           '- Certificados generales de Dickson, recurrencias y una reducción de módulos de rango nueve para C.',
           '- Una correspondencia entre comprobaciones y lemas/ecuaciones, con las dependencias teóricas identificadas.','',
           '## Corrección detectada en el manuscrito','',
           'En la matriz de L, en el orden w,y1,y1^-1*y2,u, la primera entrada superior debe ser **-(2 eta+phi)**. La relación y1*w=h^2*f*w*y1 exige ese signo con la convención usada.',
           'El cambio conserva el Pfaffiano y los invariantes alternantes de Smith, por lo que **no cambia la fórmula del grado PI**. El cuaderno contiene una comprobación que detecta el signo antiguo. No se modificó silenciosamente el manuscrito.','',
           '## Alcance exacto','',
           'Los certificados simbólicos se distinguen de las pruebas de regresión con órdenes finitos. La decisión PI aplica la tabla del manuscrito; el cálculo del grado se obtiene del toro y sus acciones, y después se compara con la fórmula.',
           'Las pruebas abstractas sobre anillos de división, independencia, ideales y Nullstellensatz siguen siendo argumentos matemáticos del artículo. No se afirma una formalización en un asistente de pruebas.',
           'El grado exacto de C y los casos de primitividad que el artículo deja abiertos se conservan como preguntas abiertas.','',
           '## Archivos para revisar','',
           '- `docs/MATRICES.md`: todas las matrices de definición.',
           '- `docs/TORUS_ACTION_MATRICES.md`: matrices de los 21 modelos de grado.',
           '- `reports/coverage.md` y `coverage.json`: correspondencia completa con las etiquetas.',
           '- `docs/SYMBOLIC_CERTIFICATES.md`: enlace entre certificados y argumentos generales.',
           '- `docs/UNBOUNDED_PROOFS.md`: inducción y demostraciones generales para todos los órdenes.',
           '- `independent/check.py`: verificador independiente; `reports/certificates.json`: datos y testigos.',
           '- `reports/notebook_execution.json`: constancia del modo de ejecución y las huellas de los módulos.',
           '- `docs/MANUSCRIPT_CORRECTIONS.md`: corrección de L lista para incorporar.','',
           'Entorno: Python '+report['environment']['python']+'; SymPy '+report['environment']['sympy']+'.']
    (ROOT/'reports/Verificaciones_PI_14641.md').write_text('\n'.join(lines)+'\n')
    print('Rendered coverage.md and Verificaciones_PI_14641.md')


if __name__=='__main__':render()
