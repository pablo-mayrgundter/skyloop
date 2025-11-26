import {Canvas} from '@react-three/fiber'
import LoopTube from '../../three/loop/LoopTube'
import Globe from '../../three/Globe'

const TubeModule = () => (
  <div className="panel">
    <h3>Tube Lab</h3>
    <p>Prototype tube geometry and anchoring points.</p>
    <div className="canvas-wrapper" style={{ height: '320px' }}>
      <Canvas camera={{ position: [0.05, 0.05, 0.15] }}>
        <ambientLight intensity={0.7} />
        <Globe />
        <LoopTube segments={128} />
      </Canvas>
    </div>
  </div>
)

export default TubeModule
