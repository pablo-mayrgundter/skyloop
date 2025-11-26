import {useFrame} from '@react-three/fiber'
import {useMemo, useRef} from 'react'
import {Group, Vector3} from 'three'
import {LOS_ANGELES, SAN_FRANCISCO} from '../../config/simulation'
import {useSimulationStore} from '../../state/simulationStore'
import {buildGreatCircleArc} from '../Globe'

type LoopGlidersProps = {
  count?: number
}

const LoopGliders = ({ count = 3 }: LoopGlidersProps) => {
  const path = useMemo(() => buildGreatCircleArc(LOS_ANGELES, SAN_FRANCISCO, 256), [])
  const gliderRefs = useRef<Group[]>([])
  const simSpeed = useSimulationStore((s) => s.simSpeed)
  const time = useSimulationStore((s) => s.time)

  useFrame((_, delta) => {
    gliderRefs.current.forEach((group, idx) => {
      const progress = ((time + delta * simSpeed) * 0.01 + idx / count) % 1
      const pathIndex = Math.floor(progress * (path.length - 1))
      const next = path[pathIndex]
      if (group && next) {
        group.position.copy(next)
        const lookAt = path[(pathIndex + 1) % path.length] ?? next.clone().add(new Vector3(0, 0.01, 0))
        group.lookAt(lookAt)
      }
    })
  })

  return (
    <group>
      {new Array(count).fill(0).map((_, idx) => (
        <group ref={(ref) => (gliderRefs.current[idx] = ref!)} key={idx}>
          <mesh position={[0, 0, 0]}>
            <capsuleGeometry args={[0.006, 0.02, 4, 8]} />
            <meshStandardMaterial color="#7ef3c2" emissive="#1b5e4a" />
          </mesh>
        </group>
      ))}
    </group>
  )
}

export default LoopGliders
