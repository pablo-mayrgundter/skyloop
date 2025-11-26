import { useRef } from 'react'
import { GridHelper } from 'three'
import { useFrame } from '@react-three/fiber'

interface GridProps {
  infiniteGrid?: boolean
  cellSize?: number
  sectionSize?: number
  sectionColor?: string
}

const Grid = ({ infiniteGrid = false, cellSize = 1, sectionSize = 10, sectionColor = '#888888' }: GridProps) => {
  const gridRef = useRef<GridHelper>(null)

  useFrame(() => {
    if (gridRef.current && infiniteGrid) {
      // For infinite grid effect, you could update the grid position based on camera
      // For now, we'll just create a large grid
    }
  })

  return <gridHelper ref={gridRef} args={[sectionSize * 10, sectionSize * 10, sectionColor, sectionColor]} />
}

export default Grid

