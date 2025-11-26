import {Canvas} from '@react-three/fiber'
import LoopTethers from '../../three/loop/LoopTethers'
import Globe from '../../three/Globe'

const TethersModule = () => (
  <div className="panel">
    <h3>Tethers Lab</h3>
    <p>Visualize tether layout and tension gradients.</p>
    <div className="canvas-wrapper" style={{ height: '320px' }}>
      <Canvas camera={{ position: [0.08, 0.08, 0.18] }}>
        <ambientLight intensity={0.5} />
        <Globe />
        <LoopTethers count={8} lengthKm={25} />
      </Canvas>
    </div>
  </div>
)

export default TethersModule
