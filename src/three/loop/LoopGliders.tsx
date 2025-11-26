import { useFrame } from '@react-three/fiber'
import { useMemo, useRef } from 'react'
import { Mesh, Vector3 } from 'three'
import { EARTH_RENDER_RADIUS, LOS_ANGELES, SAN_FRANCISCO, interpolateGreatCircle } from '../../config/simulation'
import { useSimulationStore } from '../../state/simulationStore'

const LoopGliders = () => {
  const gliderRef = useRef<Mesh>(null)
  const { time, simSpeed } = useSimulationStore((state) => ({
    time: state.time,
    simSpeed: state.simSpeed
  }))

  const path = useMemo(() => {
    const pts: Vector3[] = []
    const segments = 64
    for (let i = 0; i <= segments; i += 1) {
      pts.push(interpolateGreatCircle(LOS_ANGELES, SAN_FRANCISCO, i / segments, EARTH_RENDER_RADIUS + 0.02))
    }
    return pts
  }, [])

  useFrame(() => {
    if (!gliderRef.current) return
    const loopTime = (time * simSpeed) % path.length
    const idx = Math.floor(loopTime)
    const t = loopTime - idx
    const current = path[idx % path.length]
    const next = path[(idx + 1) % path.length]
    const position = current.clone().lerp(next, t)
    gliderRef.current.position.copy(position)
    gliderRef.current.lookAt(next)
  })

  return (
    <mesh ref={gliderRef}>
      <capsuleGeometry args={[0.01, 0.05, 8, 16]} />
      <meshStandardMaterial color="#ffdf6b" emissive="#ffbe3b" emissiveIntensity={0.8} />
    </mesh>
  )
}

export default LoopGliders
